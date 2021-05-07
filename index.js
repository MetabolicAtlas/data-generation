const fs = require('fs'), path = require('path');
const yaml = require('js-yaml');
const func = require('./func.js');

const createCsvWriter = require("csv-writer").createObjectCsvWriter;
let csvWriter = null;

let extNodeIdTracker = 1;
const humanGeneIdSet = new Set();
const externalIdDBMap = {};
const PMIDSset = new Set();
let instructions = [];
let dropIndexes = false;
let prefix = '' ;
let outputPath = '';
let outDir = './data';


const parseModelFiles = (modelDir) => {
  // find the yaml in the folder
  yamlFile = func.getFile(modelDir, /.*[.](yaml|yml)$/);
  if (!yamlFile) {
    throw new Error("yaml file not found in path ", modelDir);
  }

  const [ metadata, metabolites, reactions, genes, compartments ] = yaml.safeLoad(fs.readFileSync(yamlFile, 'utf8'));
  const metadataSection = metadata.metaData || metadata.metadata;
  const model = func.toLabelCase(metadataSection.short_name);
  // console.log("model=", model);
  const version = `V${metadataSection.version.replace(/\./g, '_')}`;
  const isHuman = metadataSection.short_name === 'Human-GEM';

  prefix = `${model}${version}`;
  outputPath = `${outDir}/${prefix}.`;

  content = { // reformat object as proper key:value objects, rename/add/remove keys
    compartmentalizedMetabolite: func.reformatCompartmentalizedMetaboliteObjets(metabolites.metabolites),
    reaction: func.reformatReactionObjets(reactions.reactions),
    gene: func.reformatGeneObjets(genes.genes),
    compartment: func.reformatCompartmentObjets(compartments.compartments),
  }

  const componentIdDict = {}; // store for each type of component the key  Id <-> element
  // use to filter out annotation/external ids for components not in the model and to add missing information
  // extracted from these annotation files such as description, etc...
  Object.keys(content).forEach((k) => {
    componentIdDict[k] = Object.fromEntries(content[k].map(e => [e[`${k}Id`], e]));
  });

  if (isHuman) {
    Object.keys(componentIdDict.gene).forEach((geneId) => {
      humanGeneIdSet.add(geneId);
    });
  }

  // subsystems are not a section in the yaml file, but are extracted from the reactions info
  content.subsystem = [];
  componentIdDict.subsystem = {};
  content.reaction.forEach((r) => {
    r.subsystems.forEach((name) => {
      const id = func.idfyString(name);
      const subsystemObject = { subsystemId: id, name }; // TODO add 'description' key
      if (!(id in componentIdDict.subsystem)) {
        content.subsystem.push(subsystemObject);
        componentIdDict.subsystem[id] = subsystemObject;
      };
    });
  });

  // ========================================================================
  // SVG mapping file
  const svgNodes = [];
  ['compartment', 'subsystem', 'custom'].forEach((component) => {
    func.createComponentSVGMapFile(component, outputPath, svgNodes, modelDir);
  });

  // write svgMaps file
  csvWriter = createCsvWriter({
    path: `${outputPath}svgMaps.csv`,
    header: svgNodes.length ? Object.keys(svgNodes[0]).map(k => Object({ id: k, title: k })) : '',
  });
  csvWriter.writeRecords(svgNodes);

  // ========================================================================
  // external IDs and annotation
  // extract EC code and PMID from YAML file
  func.createPMIDFile(PMIDSset, componentIdDict, outputPath);

  // extract information from gene annotation file
  func.extractGeneAnnotation(componentIdDict, modelDir);

  // extract description subsystem annotation file
  // TODO or remove annotation file

  // ========================================================================
  // parse External IDs files
  const externalIdNodes = [];

  ['reaction', 'metabolite', 'gene', 'subsystem'].forEach((component) => {
    extNodeIdTracker = func.createComponentExternalDbFile(externalIdNodes, externalIdDBMap,
      extNodeIdTracker, component, componentIdDict, modelDir, outputPath);
  });

  // write the externalDbs file
  csvWriter = createCsvWriter({
    path: `${outputPath}externalDbs.csv`,
    header: [{ id: 'id', title:'id' },
             { id: 'dbName', title:'dbName' },
             { id: 'externalId', title:'externalId' },
             { id: 'url', title:'url' }],
  });
  csvWriter.writeRecords(externalIdNodes);

  // ========================================================================
  // write main nodes relationships files
  // need a map to get the compartment ID from the compartment letter
  const compartmentLetterToIdMap = content.compartment.reduce((entries, c) => {
    return {
      ...entries,
      [c.letterCode]: c.compartmentId,
    };
  }, {});

  csvWriter = createCsvWriter({
    path: `${outputPath}compartmentalizedMetaboliteCompartments.csv`,
    header: [{ id: 'compartmentalizedMetaboliteId', title: 'compartmentalizedMetaboliteId' }, { id: 'compartmentId', title: 'compartmentId' }],
  });

  csvWriter.writeRecords(content.compartmentalizedMetabolite.map(
    (e) => { return { compartmentalizedMetaboliteId: e.compartmentalizedMetaboliteId, compartmentId: compartmentLetterToIdMap[e.compartment] }; }
  ));

  // ========================================================================
  // write metabolite-compartmentalizedMetabolite relationships
  // generate unique metabolite
  // keep only distinct metabolite (non-compartmentalize) and use the name to generate IDs
  let hm = {}
  const uniqueCompartmentalizedMap = {}
  content.compartmentalizedMetabolite.forEach((m) => {
    func.getUniqueCompartmentlizedMap(m, hm, uniqueCompartmentalizedMap);
  })

  const uniqueMetDict = {};
  const uniqueMetabolites = [];
  content.compartmentalizedMetabolite.forEach((m) => {
    func.getUniqueMetabolite(m, uniqueCompartmentalizedMap, uniqueMetDict, uniqueMetabolites);
  })

  // create compartmentalizedMetabolite file
  csvWriter = createCsvWriter({
    path: `${outputPath}compartmentalizedMetabolites.csv`,
    header: [{ id: 'id', title: 'id' }],
  });

  csvWriter.writeRecords(content.compartmentalizedMetabolite.map(
    (e) => { return { id: e.compartmentalizedMetaboliteId }; }
  ));

  // ========================================================================
  // extract information from metabolite annotation file
  // METABOLITES.tsv has been removed for the format, and this file is actually
  // empty in the old format

  // ========================================================================
  // CM-M relationships
  csvWriter = createCsvWriter({
    path: `${outputPath}compartmentalizedMetaboliteMetabolites.csv`,
    header: [{ id: 'compartmentalizedMetaboliteId', title: 'compartmentalizedMetaboliteId' }, { id: 'metaboliteId', title: 'metaboliteId' }],
  });

  csvWriter.writeRecords(content.compartmentalizedMetabolite.map(
    (e) => { 
      return { compartmentalizedMetaboliteId: e.compartmentalizedMetaboliteId,
               metaboliteId: uniqueCompartmentalizedMap[e.compartmentalizedMetaboliteId] }; }
  ));

  // delete compartmentlizedMetabolites, add unique metabolites
  content.metabolite = uniqueMetabolites;
  delete content.compartmentalizedMetabolite;

  // write reactants-reaction, reaction-products, reaction-genes, reaction-susbsystems relationships files
  csvWriterRR = createCsvWriter({
    path: `${outputPath}compartmentalizedMetaboliteReactions.csv`,
    header: [{ id: 'compartmentalizedMetaboliteId', title: 'compartmentalizedMetaboliteId' },
             { id: 'reactionId', title: 'reactionId' },
             { id: 'stoichiometry', title: 'stoichiometry' }],
  });
  csvWriterRP = createCsvWriter({
    path: `${outputPath}reactionCompartmentalizedMetabolites.csv`,
    header: [{ id: 'reactionId', title: 'reactionId' },
             { id: 'compartmentalizedMetaboliteId', title: 'compartmentalizedMetaboliteId' },
             { id: 'stoichiometry', title: 'stoichiometry' }],
  });
  csvWriterRG = createCsvWriter({
    path: `${outputPath}reactionGenes.csv`,
    header: [{ id: 'reactionId', title: 'reactionId' },
             { id: 'geneId', title: 'geneId' }],
  });
  csvWriterRS = createCsvWriter({
    path: `${outputPath}reactionSubsystems.csv`,
    header: [{ id: 'reactionId', title: 'reactionId' },
             { id: 'subsystemId', title: 'subsystemId' }],
  });

  const reactionReactantRecords = [];
  const reactionProductRecords = [];
  const reactionGeneRecords = [];
  const reactionSubsystemRecords = [];
  content.reaction.forEach((r) => {
    Object.entries(r.metabolites).forEach((e) => {
      const [ compartmentalizedMetaboliteId, stoichiometry ] = e;
      if (stoichiometry < 0) {
        reactionReactantRecords.push({ compartmentalizedMetaboliteId, reactionId: r.reactionId, stoichiometry: -stoichiometry });
      } else {
        reactionProductRecords.push({ reactionId: r.reactionId, compartmentalizedMetaboliteId, stoichiometry });
      }
    });
    func.getGeneIdsFromGeneRule(r.geneRule).forEach((geneId) => {
      reactionGeneRecords.push({ reactionId: r.reactionId, geneId });
    });
    r.subsystems.forEach((name) => {
      reactionSubsystemRecords.push({ reactionId: r.reactionId, subsystemId: func.idfyString(name) });
    })
  });

  csvWriterRR.writeRecords(reactionReactantRecords);
  csvWriterRP.writeRecords(reactionProductRecords);
  csvWriterRG.writeRecords(reactionGeneRecords);
  csvWriterRS.writeRecords(reactionSubsystemRecords);

  // ========================================================================
  // write nodes files
  Object.keys(content).forEach((k) => {
    const elements = content[k];
    csvWriter = createCsvWriter({
      path: `${outputPath}${k}s.csv`,
      header: [Object({ id: 'id', title: 'id' })],
    });
    csvWriter.writeRecords(elements.map(e => Object({ id: e[`${k}Id`] })));
    csvWriter = createCsvWriter({
      path: `${outputPath}${k}States.csv`,
      header: Object.keys(elements[0]).
        // ignore some keys 'metabolites', 'subsystems' are in reactions, 'compartment' is in metabolite
        filter(k => !['metabolites', 'subsystems', 'compartment'].includes(k)).
        map(k => Object({ id: k, title: k })),
    });
    // destructure object to remove the keys
    csvWriter.writeRecords(elements.map(({ subsystems, metabolites, compartment, ...e }) => e));
  });

  // TODO generate instructions more dynamically
  instructions = func.getModelCypherInstructions(prefix, dropIndexes, model, version, instructions);
};


const args = [];
try {
  for (let i = 0; i < process.argv.length; i += 1) {
    if (process.argv[i] === "--reset-db") {
      dropIndexes = true;
    } else {
      args.push(process.argv[i]);
    }
  }
  inputDir = args[2] + '/integrated-models';
} catch {
  console.log("Usage: yarn start input_dir");
  console.log("Usage: yarn start input_dir --reset-db");
  return;
}

if (!fs.existsSync('./data')){
  fs.mkdirSync('./data');
}

try {
  const intputDirFiles = fs.readdirSync(inputDir);
  for(let i = 0; i < intputDirFiles.length; i++) {
    const filePath = path.join(inputDir, intputDirFiles[i]);
    const stat = fs.lstatSync(filePath);
    if (stat.isDirectory() && intputDirFiles[i][0] != '.') {
      parseModelFiles(filePath);
    }
  }
  instructions = func.getRemainCyperInstructions(instructions);
} catch (e) {
  if (e.mark) {
    // avoid to print the whole yaml into console
    e.mark.buffer = '';
  }
  console.log(e);
  return;
}

// write cyper intructions to file
// const cyperFile = `${outDir}/import.cypher`;
fs.writeFileSync(`${outDir}/import.cypher`, instructions.join('\n'), 'utf8');

  // ========================================================================
  // write a smaller version of the hpa rna levels file, to send to the frontend
  // remove expressions of genes not in any human models parsed
if (!fs.existsSync(`${inputDir}/hpaRnaFull.json`)) {
    throw new Error("HPA rna JSON file not found");
} else {
  const hpaRnaExpressionJson = require(`${inputDir}/hpaRnaFull.json`);

  Object.keys(hpaRnaExpressionJson.levels).forEach((geneId) => {
    if (!humanGeneIdSet.has(geneId)) {
      delete hpaRnaExpressionJson.levels[geneId];
    }
  });

  const json_rna = JSON.stringify(hpaRnaExpressionJson);
  fs.writeFileSync('./data/hpaRna.json', json_rna);
}

