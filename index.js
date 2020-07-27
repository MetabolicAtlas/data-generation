
const fs = require('fs'), path = require('path');
const yaml = require('js-yaml');

const getFile = (dirPath, regexpOrString) => {
  if (!fs.existsSync(dirPath)){
    console.log("Error: no dir ", dirPath);
    return;
  }

  const files = fs.readdirSync(dirPath);
  for(let i = 0; i < files.length; i++) {
    const filePath = path.join(dirPath, files[i]);
    const stat = fs.lstatSync(filePath);
    if (!stat.isDirectory() && (regexpOrString === files[i] || (regexpOrString.test && regexpOrString.test(files[i])))) {
      return filePath;
    }
  }
};

const createCsvWriter = require("csv-writer").createObjectCsvWriter;
let csvWriter = null;

const idfyString = s => s.toLowerCase().replace(/[^a-z0-9]/g, '_').replace(/_+/g, '_').replace(/^_|_$/, ''); // for subsystems, compartments etc..
const idfyString2 = s => s.toLowerCase().replace(/[^a-z0-9]/g, '_').replace(/_+/g, '_'); // to generate compartmentalizedMetabolite ID from their name

const toLabelCase = (modelName) =>
  modelName.replace('-', ' ').split(/\s/g).map(word => `${word[0].toUpperCase()}${word.slice(1).toLowerCase()}`).join('');

const mergedObjects = data => data.reduce((acc, item) => {
  const [key, value] = Object.entries(item)[0];
    return {
    ...acc,
        [key]: value,
    };
}, {});

const getGeneIdsFromGeneRule = (geneRule) => {
  let idList = [];
  if (geneRule) {
    idList = geneRule.split(/[\s+and\s+|\s+or\s+|(+|)+|\s+]/).filter(e => e);
  }
  return idList;
}

const reformatCompartmentObjets = (data) => {
  return data.map((c) => {
    name = Object.values(c)[0];
    return { compartmentId: idfyString(name), name, letterCode: Object.keys(c)[0] };
  });
};

const reformatGeneObjets = (data) => {
  return data.map((g) => {
    id = Object.values(g[0])[0];
    return { geneId: id, name: '', alternateName: '', synonyms: '', function: '' };
  });
};

const reformatCompartmentalizedMetaboliteObjets = (data) => {
  return data.map((m) => {
    m = mergedObjects(m);
    return {
      compartmentalizedMetaboliteId: m.id,
      name: m.name,
      alternateName: '',
      synonyms: '',
      description: '',
      formula: m.formula,
      charge: m.charge,
      isCurrency: false,
      compartment: m.compartment,
    };
  });
};

const reformatReactionObjets = (data) => {
  return data.map((r) => {
    // reactionId,name,reversible,lowerBound,upperBound,geneRule,ec
    r = mergedObjects(r);
    r.metabolites = mergedObjects(r.metabolites);
    return {
      reactionId: r.id,
      name: r.name,
      metabolites: r.metabolites,
      lowerBound: r.lower_bound,
      upperBound: r.upper_bound,
      geneRule: r.gene_reaction_rule,
      ec: r.eccodes,
      subsystems: r.subsystem ? Array.isArray(r.subsystem) ? r.subsystem : [r.subsystem] : [],
    };
  } );
};

let extNodeIdTracker = 1
const humanGeneIdSet = new Set();
const externalIdDBMap = {};
const PMIDSset = new Set();
let instructions = [];
let dropIndexes = false;

const parseModelFiles = (modelDir) => {
  // find the yaml in the folder
  yamlFile = getFile(modelDir, /.*[.](yaml|yml)$/);
  if (!yamlFile) {
    console.log("Error: yaml file not found in path ", modelDir);
    return;
  }
  console.log("Files created:");

  const [ metadata, metabolites, reactions, genes, compartments ] = yaml.safeLoad(fs.readFileSync(yamlFile, 'utf8'));
  const metadataSection = metadata.metaData || metadata.metadata;
  const model = toLabelCase(metadataSection.short_name);
  const version = `V${metadataSection.version.replace(/\./g, '_')}`;
  const isHuman = metadataSection.organism === 'Homo sapiens';

  const prefix = `${model}${version}`;
  const outputPath = `./data/${prefix}.`;

  content = { // reformat object as proper key:value objects, rename/add/remove keys
    compartmentalizedMetabolite: reformatCompartmentalizedMetaboliteObjets(metabolites.metabolites),
    reaction: reformatReactionObjets(reactions.reactions),
    gene: reformatGeneObjets(genes.genes),
    compartment: reformatCompartmentObjets(compartments.compartments),
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
      const id = idfyString(name);
      const subsystemObject = { id, name }; // TODO add 'description' key
      if (!(id in componentIdDict.subsystem)) {
        content.subsystem.push(subsystemObject);
        componentIdDict.subsystem[id] = subsystemObject;
      };
    });
  });

  // ========================================================================
  // SVG mapping file
  const svgNodes = [];
  ['compartment', 'subsystem'].forEach((component) => {
    const filename = `${component}SVG.tsv`;
    const mappingFile = getFile(modelDir, filename);

    let svgRels = [];
    if (mappingFile) {
      let lines = fs.readFileSync(mappingFile, 
            { encoding: 'utf8', flag: 'r' }).split('\n').filter(Boolean);
      const filenameSet = new Set(); // check uniqness of values in the file
      for (let i = 0; i < lines.length; i++) {
        if (lines[i][0] == '#' || lines[i][0] == '@') {
          continue;
        }
        const [ componentName, mapName, mapFilename ] = lines[i].split('\t').map(e => e.trim());

        if (!content[component].map(e => e.name).includes(componentName)) {
          console.log(`Error: ${componentName} ${component} does not exist in the model`);
          exit;
        }

        if (filenameSet.has(mapFilename)) {
          console.log(`Error: map ${mapFilename} can only be linked to one ${component}`);
          exit;
        }
        filenameSet.add(mapFilename)

        if (!/^[a-z0-9_]+[.]svg$/.test(mapFilename)) {
          console.log(`Error: map ${mapFilename} (${filename}) is invalid`);
          exit;
        }
        svgNodes.push({ id: mapFilename.split('.')[0], filename: mapFilename, customName: mapName });
        svgRels.push({ [`${component}Id`]: idfyString(componentName), svgMapId: mapFilename.split('.')[0]});
      }
    } else {
      console.log(`Warning: cannot find mappingfile ${filename} in path`, modelDir);
    }

    // write the associated file
    csvWriter = createCsvWriter({
      path: `${outputPath}${component}SvgMaps.csv`,
      header: [{ id: `${component}Id`, title: `${component}Id` },
               { id: 'svgMapId', title: 'svgMapId' }],
    });
    csvWriter.writeRecords(svgRels).then(() => {
      console.log(`${component}SvgMaps.csv`);
    });
  });

  // write svgMaps file
  csvWriter = createCsvWriter({
    path: `${outputPath}svgMaps.csv`,
    header: svgNodes.length ? Object.keys(svgNodes[0]).map(k => Object({ id: k, title: k })) : '',
  });
  csvWriter.writeRecords(svgNodes).then(() => {
    console.log(`svgMaps.csv`);
  });

  // ========================================================================
  // external IDs and annotation

  // extract EC code and PMID from reaction annotation file
  const reactionAnnoFile = getFile(modelDir, /REACTIONS[.]tsv$/);
  if (!reactionAnnoFile) {
    console.log("Error: reaction annotation file not found in path", modelDir);
    return;
  }

  // TODO use one of the csv parsing lib (sync)
  lines = fs.readFileSync(reactionAnnoFile, 
            { encoding: 'utf8', flag: 'r' }).split('\n').filter(Boolean);
  const reactionPMID = [];
  const PMIDs = [];
  for (let i = 0; i < lines.length; i++) {
    if (lines[i][0] == '#' || lines[i][0] == '@') {
      continue;
    }
    const [ reactionId, ECList, PMIDList ] = lines[i].split('\t').map(e => e.trim());
    // EC are already provided by the YAML (without the 'EC:' prefix), TODO remove from annotation file?
    if (reactionId in componentIdDict.reaction && PMIDList) { //only keep the ones in the model
      PMIDList.split('; ').forEach((pubmedReferenceId) => {
        reactionPMID.push({ reactionId, pubmedReferenceId });
        if (!PMIDSset.has(pubmedReferenceId)) {
          PMIDs.push(pubmedReferenceId);
          PMIDSset.add(pubmedReferenceId);
        }
      });
    }
  }

  // create pubmedReferences file
  csvWriter = createCsvWriter({
    path: `${outputPath}pubmedReferences.csv`,
    header: [{ id: 'id', title: 'id' }],
  });
  csvWriter.writeRecords(PMIDs.map(
    (id) => { return { id }; }
  )).then(() => {
    console.log('pubmedReferences');
  });

  // write reaction pubmed reference file
  csvWriter = createCsvWriter({
    path: `${outputPath}reactionPubmedReferences.csv`,
    header: [{ id: 'reactionId', title: 'reactionId' },
             { id: 'pubmedReferenceId', title: 'pubmedReferenceId' }],
  });
  csvWriter.writeRecords(reactionPMID).then(() => {
    console.log('reactionPubmedReferences');
  });

  // extract information from gene annotation file

  const geneAnnoFile = getFile(modelDir, /GENES[.]tsv$/);
  if (!geneAnnoFile) {
    console.log("Error: gene annotation file not found in path", modelDir);
    return;
  }

  // TODO use one of the csv parsing lib (sync)
  lines = fs.readFileSync(geneAnnoFile, 
            { encoding: 'utf8', flag: 'r' }).split('\n').filter(Boolean);
  for (let i = 0; i < lines.length; i++) {
    if (lines[i][0] == '#' || lines[i][0] == '@') {
      continue;
    }
    const [ geneId, name, alternateName, synonyms, thefunction, ec, catalytic_activity ] = lines[i].split('\t').map(e => e.trim());
    if (geneId in componentIdDict.gene) { //only keep the ones in the model
      const gene = componentIdDict.gene[geneId];
      Object.assign(gene, { name, alternateName, synonyms, function: thefunction }); // other props are not in the db design, TODO remove them?
    }
  }

  // extract description subsystem annotation file
  // TODO or remove annotation file

  // ========================================================================
  // parse External IDs files
  const externalIdNodes = [];

  ['reaction', 'metabolite', 'gene', 'subsystem'].forEach((component) => {
    const externalIdDBComponentRel = [];
    const filename = `${component.toUpperCase()}S_EID.tsv`;
    const extIDFile = getFile(modelDir, filename);
    const IdSetKey = component === 'metabolite' ? 'compartmentalizedMetabolite' : component;

    if (extIDFile) {
      // TODO use one of the csv parsing lib (sync)
      lines = fs.readFileSync(extIDFile, 
                { encoding: 'utf8', flag: 'r' }).split('\n').filter(Boolean);

      for (let i = 0; i < lines.length; i++) {
        if (lines[i][0] == '#' || lines[i][0] == '@') {
          continue;
        }
        const [ id, dbName, externalId, url ] = lines[i].split('\t').map(e => e.trim());
        if (!(id in componentIdDict[IdSetKey])) { //only keep the ones in the model
          continue;
        }

        const externalDbEntryKey = `${dbName}${externalId}${url}`; // diff url leads to new nodes!

        let node = null;
        if (externalDbEntryKey in externalIdDBMap) {
          node = externalIdDBMap[externalDbEntryKey]; // reuse the node and id
        } else {
          node = { id: extNodeIdTracker, dbName, externalId, url };
          externalIdDBMap[externalDbEntryKey] = node;
          extNodeIdTracker += 1;

          // save the node for externalDBs.csv
          externalIdNodes.push(node);
        }

        // save the relationships between the node and the current component ID (reaction, gene, etc)
        externalIdDBComponentRel.push({ id, externalDbId: node.id }); // e.g. geneId, externalDbId
      }
    } else {
      console.log(`Warning: cannot find external ID file ${filename} in path`, modelDir);
    }

    // write the associated file
    csvWriter = createCsvWriter({
      path: `${outputPath}${component}ExternalDbs.csv`,
      header: [{ id: `${component}Id`, title: `${component}Id` },
               { id: 'externalDbId', title: 'externalDbId' }],
    });
    csvWriter.writeRecords(externalIdDBComponentRel.map(
      (e) => { return { [`${component}Id`]: e.id, externalDbId: e.externalDbId }; }
    )).then(() => {
      console.log(`${component}ExternalDbs.csv`);
    });
  });

  if (externalIdNodes.length !== 0) {
    // write the externalDbs file
    csvWriter = createCsvWriter({
      path: `${outputPath}externalDbs.csv`,
      header: Object.keys(externalIdNodes[0]).map(k => Object({ id: k, title: k })),
    });
    csvWriter.writeRecords(externalIdNodes).then(() => {
      console.log('externalDbs');
    });
  }

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
  )).then(() => {
    console.log('compartmentalizedMetaboliteCompartments');
  });

  // ========================================================================
  // write metabolite-compartmentalizedMetabolite relationships
  // generate unique metabolite
  // keep only distinct metabolite (non-compartmentalize) and use the name to generate IDs
  let hm = {}
  const uniqueCompartmentalizedMap = {}
  content.compartmentalizedMetabolite.forEach((m) => {
    const newID = idfyString2(m.name);
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
  )).then(() => {
    console.log('compartmentalizedMetabolites');
  });

  // ========================================================================
  // extract information from metabolite annotation file

  const metaboliteAnnoFile = getFile(modelDir, /METABOLITES[.]tsv$/);
  if (!metaboliteAnnoFile) {
    console.log("Error: metabolite annotation file not found in path", modelDir);
    return;
  }

  // TODO use one of the csv parsing lib (sync)
  lines = fs.readFileSync(metaboliteAnnoFile, 
            { encoding: 'utf8', flag: 'r' }).split('\n').filter(Boolean);
  for (let i = 0; i < lines.length; i++) {
    if (lines[i][0] == '#' || lines[i][0] == '@') {
      continue;
    }
    const [ metaboliteId, alternateName, synonyms, description, mass, inchi ] = lines[i].split('\t').map(e => e.trim());
    if (metaboliteId in componentIdDict.compartmentalizedMetabolite) { //only keep the ones in the model
      // find the unique met associated
      const umet = uniqueMetDict[componentIdDict.compartmentalizedMetabolite[metaboliteId].name];
      Object.assign(umet, { alternateName, synonyms, description }); // other props are not in the db design, TODO remove them?
    }
  }

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
  )).then(() => {
    console.log('compartmentalizedMetaboliteMetabolites');
  });

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
    getGeneIdsFromGeneRule(r.geneRule).forEach((geneId) => {
      reactionGeneRecords.push({ reactionId: r.reactionId, geneId });
    });
    r.subsystems.forEach((name) => {
      reactionSubsystemRecords.push({ reactionId: r.reactionId, subsystemId: idfyString(name) });
    })
  });

  csvWriterRR.writeRecords(reactionReactantRecords).then(() => {
    console.log('compartmentalizedMetaboliteReactions');
  });
  csvWriterRP.writeRecords(reactionProductRecords).then(() => {
    console.log('reactionCompartmentalizedMetabolites');
  });
  csvWriterRG.writeRecords(reactionGeneRecords).then(() => {
    console.log('reactionGenes');
  });
  csvWriterRS.writeRecords(reactionSubsystemRecords).then(() => {
    console.log('reactionSubsystems');
  });

  // ========================================================================
  // write nodes files
  Object.keys(content).forEach((k) => {
    const elements = content[k];
    csvWriter = createCsvWriter({
      path: `${outputPath}${k}s.csv`,
      header: [Object({ id: 'id', title: 'id' })],
    });
    csvWriter.writeRecords(elements.map(e => Object({ id: e[`${k}Id`] }))).then(() => {
      console.log(`${k}s`);
    });
    csvWriter = createCsvWriter({
      path: `${outputPath}${k}States.csv`,
      header: Object.keys(elements[0]).
        // ignore some keys 'metabolites', 'subsystems' are in reactions, 'compartment' is in metabolite
        filter(k => !['metabolites', 'subsystems', 'compartment'].includes(k)).
        map(k => Object({ id: k, title: k })),
    });
    // destructure object to remove the keys
    csvWriter.writeRecords(elements.map(({ subsystems, metabolites, compartment, ...e }) => e)).then(() => {
      console.log(`${k}States file generated.`);
    });
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
        'CALL apoc.schema.assert({},{},true) YIELD label, key RETURN *;',
        'CALL  db.index.fulltext.drop(\"fulltext\");\n',
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
MATCH (n:Metabolite {id: csvLine.metaboliteId})
CREATE (ns:MetaboliteState {name:csvLine.name,alternateName:csvLine.alternateName,synonyms:csvLine.synonyms,description:csvLine.description,formula:csvLine.formula,charge:toInteger(csvLine.charge),isCurrency:toBoolean(csvLine.isCurrency)})
CREATE (n)-[:${version}]->(ns);

LOAD CSV WITH HEADERS FROM "file:///${prefix}.compartmentalizedMetabolites.csv" AS csvLine
CREATE (n:CompartmentalizedMetabolite:${model} {id:csvLine.id});

LOAD CSV WITH HEADERS FROM "file:///${prefix}.compartments.csv" AS csvLine
CREATE (n:Compartment:${model} {id:csvLine.id});
LOAD CSV WITH HEADERS FROM "file:///${prefix}.compartmentStates.csv" AS csvLine
MATCH (n:Compartment {id: csvLine.compartmentId})
CREATE (ns:CompartmentState {name:csvLine.name,letterCode:csvLine.letterCode})
CREATE (n)-[:${version}]->(ns);

LOAD CSV WITH HEADERS FROM "file:///${prefix}.reactions.csv" AS csvLine
CREATE (n:Reaction:${model} {id:csvLine.id});
LOAD CSV WITH HEADERS FROM "file:///${prefix}.reactionStates.csv" AS csvLine
MATCH (n:Reaction {id: csvLine.reactionId})
CREATE (ns:ReactionState {name:csvLine.name,reversible:toBoolean(csvLine.reversible),lowerBound:toInteger(csvLine.lowerBound),upperBound:toInteger(csvLine.upperBound),geneRule:csvLine.geneRule,ec:csvLine.ec})
CREATE (n)-[:${version}]->(ns);

LOAD CSV WITH HEADERS FROM "file:///${prefix}.genes.csv" AS csvLine
CREATE (n:Gene:${model} {id:csvLine.id});
LOAD CSV WITH HEADERS FROM "file:///${prefix}.geneStates.csv" AS csvLine
MATCH (n:Gene {id: csvLine.geneId})
CREATE (ns:GeneState {name:csvLine.name,alternateName:csvLine.alternateName,synonyms:csvLine.synonyms,function:csvLine.function})
CREATE (n)-[:${version}]->(ns);

LOAD CSV WITH HEADERS FROM "file:///${prefix}.subsystems.csv" AS csvLine
CREATE (n:Subsystem:${model} {id:csvLine.id});
LOAD CSV WITH HEADERS FROM "file:///${prefix}.subsystemStates.csv" AS csvLine
MATCH (n:Subsystem {id: csvLine.subsystemId})
CREATE (ns:SubsystemState {name:csvLine.name})
CREATE (n)-[:${version}]->(ns);

LOAD CSV WITH HEADERS FROM "file:///${prefix}.svgMaps.csv" AS csvLine
CREATE (n:SvgMap:${model} {id:csvLine.id,filename:csvLine.filename,customName:csvLine.customName});

LOAD CSV WITH HEADERS FROM "file:///${prefix}.externalDbs.csv" AS csvLine
CREATE (n:ExternalDb {id:csvLine.id,dbName:csvLine.dbName,externalId:csvLine.externalId,url:csvLine.url});

LOAD CSV WITH HEADERS FROM "file:///${prefix}.pubmedReferences.csv" AS csvLine
CREATE (n:PubmedReference {id:csvLine.id,pubmedId:csvLine.pubmedId});

LOAD CSV WITH HEADERS FROM "file:///${prefix}.compartmentalizedMetaboliteMetabolites.csv" AS csvLine
MATCH (n1:CompartmentalizedMetabolite {id: csvLine.compartmentalizedMetaboliteId}),(n2:Metabolite {id: csvLine.metaboliteId})
CREATE (n1)-[:${version}]->(n2);

LOAD CSV WITH HEADERS FROM "file:///${prefix}.compartmentalizedMetaboliteCompartments.csv" AS csvLine
MATCH (n1:CompartmentalizedMetabolite {id: csvLine.compartmentalizedMetaboliteId}),(n2:Compartment {id: csvLine.compartmentId})
CREATE (n1)-[:${version}]->(n2);

LOAD CSV WITH HEADERS FROM "file:///${prefix}.compartmentalizedMetaboliteReactions.csv" AS csvLine
MATCH (n1:CompartmentalizedMetabolite {id: csvLine.compartmentalizedMetaboliteId}),(n2:Reaction {id: csvLine.reactionId})
CREATE (n1)-[:${version} {stoichiometry:toFloat(csvLine.stoichiometry)}]->(n2);

LOAD CSV WITH HEADERS FROM "file:///${prefix}.reactionCompartmentalizedMetabolites.csv" AS csvLine
MATCH (n1:Reaction {id: csvLine.reactionId}),(n2:CompartmentalizedMetabolite {id: csvLine.compartmentalizedMetaboliteId})
CREATE (n1)-[:${version} {stoichiometry:toFloat(csvLine.stoichiometry)}]->(n2);

LOAD CSV WITH HEADERS FROM "file:///${prefix}.reactionGenes.csv" AS csvLine
MATCH (n1:Reaction {id: csvLine.reactionId}),(n2:Gene {id: csvLine.geneId})
CREATE (n1)-[:${version}]->(n2);

LOAD CSV WITH HEADERS FROM "file:///${prefix}.reactionSubsystems.csv" AS csvLine
MATCH (n1:Reaction {id: csvLine.reactionId}),(n2:Subsystem {id: csvLine.subsystemId})
CREATE (n1)-[:${version}]->(n2);

LOAD CSV WITH HEADERS FROM "file:///${prefix}.reactionPubmedReferences.csv" AS csvLine
MATCH (n1:Reaction {id: csvLine.reactionId}),(n2:PubmedReference {id: csvLine.pubmedReferenceId})
CREATE (n1)-[:${version}]->(n2);

LOAD CSV WITH HEADERS FROM "file:///${prefix}.compartmentSvgMaps.csv" AS csvLine
MATCH (n1:Compartment {id: csvLine.compartmentId}),(n2:SvgMap {id: csvLine.svgMapId})
CREATE (n1)-[:${version}]->(n2);

LOAD CSV WITH HEADERS FROM "file:///${prefix}.subsystemSvgMaps.csv" AS csvLine
MATCH (n1:Subsystem {id: csvLine.subsystemId}),(n2:SvgMap {id: csvLine.svgMapId})
CREATE (n1)-[:${version}]->(n2);

LOAD CSV WITH HEADERS FROM "file:///${prefix}.metaboliteExternalDbs.csv" AS csvLine
MATCH (n1:CompartmentalizedMetabolite {id: csvLine.metaboliteId}),(n2:ExternalDb {id: csvLine.externalDbId})
CREATE (n1)-[:${version}]->(n2);

LOAD CSV WITH HEADERS FROM "file:///${prefix}.subsystemExternalDbs.csv" AS csvLine
MATCH (n1:Subsystem {id: csvLine.subsystemId}),(n2:ExternalDb {id: csvLine.externalDbId})
CREATE (n1)-[:${version}]->(n2);

LOAD CSV WITH HEADERS FROM "file:///${prefix}.reactionExternalDbs.csv" AS csvLine
MATCH (n1:Reaction {id: csvLine.reactionId}),(n2:ExternalDb {id: csvLine.externalDbId})
CREATE (n1)-[:${version}]->(n2);

LOAD CSV WITH HEADERS FROM "file:///${prefix}.geneExternalDbs.csv" AS csvLine
MATCH (n1:Gene {id: csvLine.geneId}),(n2:ExternalDb {id: csvLine.externalDbId})
CREATE (n1)-[:${version}]->(n2);

`
  cypherInstructions.split('\n').forEach(i => {
    instructions.push(i);
  });
};


const args = [];
try {
  for (let i = 0; i < process.argv.length; i += 1) {
    if (process.argv[i] === "--resset-db") {
      dropIndexes = true;
    } else {
      args.push(process.argv[i]);
    }
  }
  inputDir = args[2];
} catch {
  console.log("Usage: yarn start input_dir");
  console.log("Usage: yarn start input_dir --drop-indexes");
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
    ["CompartmentState", "Compartment", "MetaboliteState", "Metabolite", "CompartmentalizedMetabolite", "SubsystemState", "Subsystem", "ReactionState", "Reaction", "GeneState", "Gene", "PubMedReference"],
    ["id", "name", "letterCode", "alternateName", "synonyms", "description", "formula", "function", "pubMedID"])
  `.split('\n').forEach(i => {
    instructions.push(i);
  });
} catch (e) {
  console.log(e);
  return;
}

  // write cyper intructions to file
fs.writeFileSync('./data/import.cypher', instructions.join('\n'), 'utf8');

  // ========================================================================
  // write a smaller version of the hpa rna levels file, to send to the frontend
  // remove expressions of genes not in any human models parsed
if (!fs.existsSync(`${inputDir}/hpaRnaFull.json`)) {
    console.log("Error: HPA rna JSON file not found");
    return;
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

