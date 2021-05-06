const fs = require('fs'), path = require('path');
const yaml = require('js-yaml');
const func = require('./func.js');
const { dbnameDict } = require('./var');

const createCsvWriter = require("csv-writer").createObjectCsvWriter;
let csvWriter = null;

let extNodeIdTracker = 1
const humanGeneIdSet = new Set();
const externalIdDBMap = {};
const PMIDSset = new Set();
let instructions = [];
let dropIndexes = false;

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

  const prefix = `${model}${version}`;
  const outputPath = `./data/${prefix}.`;

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
    const newID = func.idfyString2(m.name);
    if (!(newID in hm)) {
      hm[newID] = m.name;
      uniqueCompartmentalizedMap[m.compartmentalizedMetaboliteId] = newID;
    } else {
      if (hm[newID] !== m.name) {
        // console.log('Error duplicated ID:' + newID + '(' + m.name + ') collision with ' + hm[newID]);
        uniqueCompartmentalizedMap[m.compartmentalizedMetaboliteId] = newID + '_';
      } else {
        uniqueCompartmentalizedMap[m.compartmentalizedMetaboliteId] = newID;
      }
    }
  })

  const uniqueMetDict = {};
  const uniqueMetabolites = [];
  content.compartmentalizedMetabolite.forEach((m) => {
    const newID = uniqueCompartmentalizedMap[m.compartmentalizedMetaboliteId];
    if (!(m.name in uniqueMetDict)) {
      const uMet = {
        metaboliteId: newID,
        name: m.name,
        alternateName: m.alternateName,
        synonyms: m.synonyms,
        description: m.description,
        formula: m.formula,
        charge: m.charge,
        isCurrency: m.isCurrency,
      };
      uniqueMetabolites.push(uMet);
      uniqueMetDict[uMet.name] = uMet;
    }
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
  if (instructions.length === 0) {
    instructions = [
      'MATCH (n) DETACH DELETE n;\n',
    ];

    if (dropIndexes) {
      instructions = instructions.concat([
        'DROP INDEX ON :Metabolite(id);',
        'DROP INDEX ON :CompartmentalizedMetabolite(id);',
        'DROP INDEX ON :Compartment(id);',
        'DROP INDEX ON :Reaction(id);',
        'DROP INDEX ON :Gene(id);',
        'DROP INDEX ON :Subsystem(id);',
        'DROP INDEX ON :SvgMap(id);',
        'DROP INDEX ON :ExternalDb(id);',
        'DROP INDEX ON :PubmedReference(id);\n',
        'CALL db.index.fulltext.drop(\"fulltext\");\n',
      ]);
    }

    instructions = instructions.concat([
      'CREATE INDEX FOR (n:Metabolite) ON (n.id);',
      'CREATE INDEX FOR (n:CompartmentalizedMetabolite) ON (n.id);',
      'CREATE INDEX FOR (n:Compartment) ON (n.id);',
      'CREATE INDEX FOR (n:Reaction) ON (n.id);',
      'CREATE INDEX FOR (n:Gene) ON (n.id);',
      'CREATE INDEX FOR (n:Subsystem) ON (n.id);',
      'CREATE INDEX FOR (n:SvgMap) ON (n.id);',
      'CREATE INDEX FOR (n:ExternalDb) ON (n.id);',
      'CREATE INDEX FOR (n:PubmedReference) ON (n.id);\n',
    ]);
  }

  const cypherInstructions = `
LOAD CSV WITH HEADERS FROM "file:///${prefix}.metabolites.csv" AS csvLine
CREATE (n:Metabolite:${model} {id:csvLine.id});
LOAD CSV WITH HEADERS FROM "file:///${prefix}.metaboliteStates.csv" AS csvLine
MATCH (n:Metabolite:${model} {id: csvLine.metaboliteId})
CREATE (ns:MetaboliteState:${model} {name:csvLine.name,alternateName:csvLine.alternateName,synonyms:csvLine.synonyms,description:csvLine.description,formula:csvLine.formula,charge:toInteger(csvLine.charge),isCurrency:toBoolean(csvLine.isCurrency)})
CREATE (n)-[:${version}]->(ns);

LOAD CSV WITH HEADERS FROM "file:///${prefix}.compartmentalizedMetabolites.csv" AS csvLine
CREATE (n:CompartmentalizedMetabolite:${model} {id:csvLine.id});

LOAD CSV WITH HEADERS FROM "file:///${prefix}.compartments.csv" AS csvLine
CREATE (n:Compartment:${model} {id:csvLine.id});
LOAD CSV WITH HEADERS FROM "file:///${prefix}.compartmentStates.csv" AS csvLine
MATCH (n:Compartment:${model} {id: csvLine.compartmentId})
CREATE (ns:CompartmentState:${model} {name:csvLine.name,letterCode:csvLine.letterCode})
CREATE (n)-[:${version}]->(ns);

LOAD CSV WITH HEADERS FROM "file:///${prefix}.reactions.csv" AS csvLine
CREATE (n:Reaction:${model} {id:csvLine.id});
LOAD CSV WITH HEADERS FROM "file:///${prefix}.reactionStates.csv" AS csvLine
MATCH (n:Reaction:${model} {id: csvLine.reactionId})
CREATE (ns:ReactionState:${model} {name:csvLine.name,reversible:toBoolean(csvLine.reversible),lowerBound:toInteger(csvLine.lowerBound),upperBound:toInteger(csvLine.upperBound),geneRule:csvLine.geneRule,ec:csvLine.ec})
CREATE (n)-[:${version}]->(ns);

LOAD CSV WITH HEADERS FROM "file:///${prefix}.genes.csv" AS csvLine
CREATE (n:Gene:${model} {id:csvLine.id});
LOAD CSV WITH HEADERS FROM "file:///${prefix}.geneStates.csv" AS csvLine
MATCH (n:Gene:${model} {id: csvLine.geneId})
CREATE (ns:GeneState:${model} {name:csvLine.name,alternateName:csvLine.alternateName,synonyms:csvLine.synonyms,function:csvLine.function})
CREATE (n)-[:${version}]->(ns);

LOAD CSV WITH HEADERS FROM "file:///${prefix}.subsystems.csv" AS csvLine
CREATE (n:Subsystem:${model} {id:csvLine.id});
LOAD CSV WITH HEADERS FROM "file:///${prefix}.subsystemStates.csv" AS csvLine
MATCH (n:Subsystem:${model} {id: csvLine.subsystemId})
CREATE (ns:SubsystemState:${model} {name:csvLine.name})
CREATE (n)-[:${version}]->(ns);

LOAD CSV WITH HEADERS FROM "file:///${prefix}.svgMaps.csv" AS csvLine
CREATE (n:SvgMap:${model} {id:csvLine.id,filename:csvLine.filename,customName:csvLine.customName});

LOAD CSV WITH HEADERS FROM "file:///${prefix}.externalDbs.csv" AS csvLine
CREATE (n:ExternalDb {id:csvLine.id,dbName:csvLine.dbName,externalId:csvLine.externalId,url:csvLine.url});

LOAD CSV WITH HEADERS FROM "file:///${prefix}.pubmedReferences.csv" AS csvLine
CREATE (n:PubmedReference {id:csvLine.id,pubmedId:csvLine.pubmedId});

LOAD CSV WITH HEADERS FROM "file:///${prefix}.compartmentalizedMetaboliteMetabolites.csv" AS csvLine
MATCH (n1:CompartmentalizedMetabolite:${model} {id: csvLine.compartmentalizedMetaboliteId}),(n2:Metabolite:${model} {id: csvLine.metaboliteId})
CREATE (n1)-[:${version}]->(n2);

LOAD CSV WITH HEADERS FROM "file:///${prefix}.compartmentalizedMetaboliteCompartments.csv" AS csvLine
MATCH (n1:CompartmentalizedMetabolite:${model} {id: csvLine.compartmentalizedMetaboliteId}),(n2:Compartment:${model} {id: csvLine.compartmentId})
CREATE (n1)-[:${version}]->(n2);

LOAD CSV WITH HEADERS FROM "file:///${prefix}.compartmentalizedMetaboliteReactions.csv" AS csvLine
MATCH (n1:CompartmentalizedMetabolite:${model} {id: csvLine.compartmentalizedMetaboliteId}),(n2:Reaction:${model} {id: csvLine.reactionId})
CREATE (n1)-[:${version} {stoichiometry:toFloat(csvLine.stoichiometry)}]->(n2);

LOAD CSV WITH HEADERS FROM "file:///${prefix}.reactionCompartmentalizedMetabolites.csv" AS csvLine
MATCH (n1:Reaction:${model} {id: csvLine.reactionId}),(n2:CompartmentalizedMetabolite:${model} {id: csvLine.compartmentalizedMetaboliteId})
CREATE (n1)-[:${version} {stoichiometry:toFloat(csvLine.stoichiometry)}]->(n2);

LOAD CSV WITH HEADERS FROM "file:///${prefix}.reactionGenes.csv" AS csvLine
MATCH (n1:Reaction:${model} {id: csvLine.reactionId}),(n2:Gene:${model} {id: csvLine.geneId})
CREATE (n1)-[:${version}]->(n2);

LOAD CSV WITH HEADERS FROM "file:///${prefix}.reactionSubsystems.csv" AS csvLine
MATCH (n1:Reaction:${model} {id: csvLine.reactionId}),(n2:Subsystem:${model} {id: csvLine.subsystemId})
CREATE (n1)-[:${version}]->(n2);

LOAD CSV WITH HEADERS FROM "file:///${prefix}.reactionPubmedReferences.csv" AS csvLine
MATCH (n1:Reaction:${model} {id: csvLine.reactionId}),(n2:PubmedReference {id: csvLine.pubmedReferenceId})
CREATE (n1)-[:${version}]->(n2);

LOAD CSV WITH HEADERS FROM "file:///${prefix}.compartmentSvgMaps.csv" AS csvLine
MATCH (n1:Compartment:${model} {id: csvLine.compartmentId}),(n2:SvgMap:${model} {id: csvLine.svgMapId})
CREATE (n1)-[:${version}]->(n2);

LOAD CSV WITH HEADERS FROM "file:///${prefix}.subsystemSvgMaps.csv" AS csvLine
MATCH (n1:Subsystem:${model} {id: csvLine.subsystemId}),(n2:SvgMap:${model} {id: csvLine.svgMapId})
CREATE (n1)-[:${version}]->(n2);

LOAD CSV WITH HEADERS FROM "file:///${prefix}.compartmentalizedMetaboliteExternalDbs.csv" AS csvLine
MATCH (n1:CompartmentalizedMetabolite:${model} {id: csvLine.compartmentalizedMetaboliteId}),(n2:ExternalDb {id: csvLine.externalDbId})
CREATE (n1)-[:${version}]->(n2);

LOAD CSV WITH HEADERS FROM "file:///${prefix}.subsystemExternalDbs.csv" AS csvLine
MATCH (n1:Subsystem:${model} {id: csvLine.subsystemId}),(n2:ExternalDb {id: csvLine.externalDbId})
CREATE (n1)-[:${version}]->(n2);

LOAD CSV WITH HEADERS FROM "file:///${prefix}.reactionExternalDbs.csv" AS csvLine
MATCH (n1:Reaction:${model} {id: csvLine.reactionId}),(n2:ExternalDb {id: csvLine.externalDbId})
CREATE (n1)-[:${version}]->(n2);

LOAD CSV WITH HEADERS FROM "file:///${prefix}.geneExternalDbs.csv" AS csvLine
MATCH (n1:Gene:${model} {id: csvLine.geneId}),(n2:ExternalDb {id: csvLine.externalDbId})
CREATE (n1)-[:${version}]->(n2);

`
  cypherInstructions.split('\n').forEach(i => {
    instructions.push(i);
  });
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
  `CALL db.index.fulltext.createNodeIndex(
    "fulltext",
    ["CompartmentState", "Compartment", "MetaboliteState", "Metabolite", "CompartmentalizedMetabolite", "SubsystemState", "Subsystem", "ReactionState", "Reaction", "GeneState", "Gene", "PubmedReference"],
    ["id", "name", "letterCode", "alternateName", "synonyms", "description", "formula", "function", "pubMedID", "ec"]);
  `.split('\n').forEach(i => {
    instructions.push(i);
  });
} catch (e) {
  if (e.mark) {
    // avoid to print the whole yaml into console
    e.mark.buffer = '';
  }
  console.log(e);
  return;
}

  // write cyper intructions to file
fs.writeFileSync('./data/import.cypher', instructions.join('\n'), 'utf8');

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

