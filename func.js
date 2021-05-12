const fs = require('fs'), path = require('path');
const createCsvWriter = require("csv-writer").createObjectCsvWriter;
const { dbnameDict } = require('./var.js');
const utils = require('./utils.js');


const createComponentSVGMapFile = (component, outputPath, svgNodes, modelDir) => {
  const filename = `${component}SVG.tsv.plain`;
  const mappingFile = utils.getFile(modelDir, filename);
  const isCustom = component === 'custom';

  let svgRels = [];
  if (mappingFile) {
    let lines = fs.readFileSync(mappingFile, 
          { encoding: 'utf8', flag: 'r' }).split('\n').filter(Boolean);
    const filenameSet = new Set(); // check uniqness of values in the file
    for (let i = 0; i < lines.length; i++) {
      if (lines[i][0] == '#' || lines[i][0] == '@') {
        continue;
      }

      let componentName, mapName, mapFilename;

      const columns = lines[i].split('\t').map(e => e.trim());
      if (isCustom) {
        [ mapName, mapFilename ] = columns;
      } else {
        [ componentName, mapName, mapFilename ] = columns;
      }

      if (componentName && !content[component].map(e => e.name).includes(componentName)) {
        throw new Error(`${component} "${componentName}" does not exist in the model "${metadataSection.short_name}"`);
      }

      if (filenameSet.has(mapFilename)) {
        throw new Error(`map ${mapFilename} can only be linked to one ${component} in the model "${metadataSection.short_name}"`);
      }
      filenameSet.add(mapFilename)

      if (!/^[a-z0-9_]+[.]svg$/.test(mapFilename)) {
        throw new Error(`map "${mapFilename}" referenced by ${metadataSection.short_name}/${filename} is invalid`);
      }
      svgNodes.push({ id: mapFilename.split('.')[0], filename: mapFilename, customName: mapName });

      if (componentName) {
        svgRels.push({
          [`${component}Id`]: utils.idfyString(componentName),
          svgMapId: mapFilename.split('.')[0],
        });
      }
    }
  } else {
    console.log(`Warning: cannot find mappingfile ${filename} in path`, modelDir);
  }

  // write the associated file
  const csvWriter = createCsvWriter({
    path: `${outputPath}${component}SvgMaps.csv`,
    header: [{ id: `${component}Id`, title: `${component}Id` },
              { id: 'svgMapId', title: 'svgMapId' }],
  });
  csvWriter.writeRecords(svgRels);
};

const createPMIDFile = (PMIDSset, componentIdDict, outputPath) => {
  const reactionPMID = [];
  const PMIDs = [];
  for (const reactionId in componentIdDict.reaction) {
    if (reactionId.match('^HMR_')) {
      const ECList = componentIdDict.reaction[reactionId].ec;
      const PMIDList = componentIdDict.reaction[reactionId].references;
      if (PMIDList) {
        PMIDList.split(';').forEach((pubmedReferenceId) => {
          if (pubmedReferenceId.match('^PMID')) {
            pubmedReferenceId = pubmedReferenceId.replace(/PMID:*/g, '');
            // console.log(pubmedReferenceId);
            reactionPMID.push({ reactionId, pubmedReferenceId });
            if (!PMIDSset.has(pubmedReferenceId)) {
              PMIDs.push(pubmedReferenceId);
              PMIDSset.add(pubmedReferenceId);
            }
          }
        });
      }
    }
  }
  // create pubmedReferences file
  csvWriter = createCsvWriter({
    path: `${outputPath}pubmedReferences.csv`,
    header: [{ id: 'id', title: 'id' }],
  });
  csvWriter.writeRecords(PMIDs.map(
    (id) => { return { id }; }
  ));

  // write reaction pubmed reference file
  csvWriter = createCsvWriter({
    path: `${outputPath}reactionPubmedReferences.csv`,
    header: [{ id: 'reactionId', title: 'reactionId' },
             { id: 'pubmedReferenceId', title: 'pubmedReferenceId' }],
  });
  csvWriter.writeRecords(reactionPMID);

};

const extractGeneAnnotation = (componentIdDict, modelDir) => {
  const geneAnnoFile = utils.getFile(modelDir, /genes-new[.]tsv$/);
  if (!geneAnnoFile) {
    console.log("Warning: cannot find gene annotation file genes-new.tsv in path", modelDir);
  } else {
    // TODO use one of the csv parsing lib (sync)
    lines = fs.readFileSync(geneAnnoFile, 
              { encoding: 'utf8', flag: 'r' }).split('\n').filter(Boolean);
    for (let i = 0; i < lines.length; i++) {
      if (lines[i][0] == '#' || lines[i][0] == '@') {
        continue;
      }
      // thefunction, ec and catalytic_activity are not defined in the new TSV
      // format and thus set the default value as empty 
      const thefunction = '';
      const ec = '';
      const catalytic_activity = '';
      const [ geneId, geneENSTID, geneENSPID, geneUniProtID, name, geneEntrezID, alternateName, synonyms] = lines[i].split('\t').map(e => utils.trim(e, '"'));
      if (geneId in componentIdDict.gene) { //only keep the ones in the model
        const gene = componentIdDict.gene[geneId];
        Object.assign(gene, { name, alternateName, synonyms, function: thefunction }); // other props are not in the db design, TODO remove them?
      }
    }
  }
}

const createComponentExternalDbFile = (externalIdNodes, externalIdDBMap, extNodeIdTracker, component, componentIdDict, modelDir, outputPath) => {
  const externalIdDBComponentRel = [];
  const filename = `${component}s-new.tsv`;
  const extIDFile = utils.getFile(modelDir, filename);
  const fcomponent = component === 'metabolite' ? 'compartmentalizedMetabolite' : component;

  if (extIDFile) {
    // TODO use one of the csv parsing lib (sync)
    lines = fs.readFileSync(extIDFile, 
              { encoding: 'utf8', flag: 'r' }).split('\n').filter(Boolean);

    var headerArr = [];
    var contentArr = [];
    for (let i = 0; i < lines.length; i++) {
      if (lines[i][0] == '#') {
        continue;
      } else if (i == 0) { /*read the header line*/
        headerArr = lines[i].split('\t').map(e => e.trim());
        continue;
      } else {
        contentArr = lines[i].split('\t').map(e => utils.trim(e, '"'));
      }

      const id = contentArr[0];
      if (!(id in componentIdDict[fcomponent])) { //only keep the ones in the model
        console.log('Warning: id ' + id + ' not in '  + ' componentIdDict[' + fcomponent+']');
        continue;
      }

      if (fcomponent == "gene"){ /*add two more items Ensembl and Protein Atlas which is not included in the new format*/
        headerArr.push('geneEnsemblID');
        headerArr.push('geneProteinAtlasID');
        contentArr.push(id); /*For Ensembl, externalId is equal to id*/
        contentArr.push(id); /*For Protein Atlas, externalId is equal to id*/  
      }
      const numItem = contentArr.length;

      for (let j = 1; j < numItem; j++) {
        const header = headerArr[j];
        const regexGene = "gene.*ID$";
        const regexRxn = "rxn.*ID$";
        const regexMet = "met.*ID$";
        if ((fcomponent == 'gene' && header.match(regexGene) == null) ||
            (fcomponent == 'reaction' && header.match(regexRxn) == null) ||
            (fcomponent == 'compartmentalizedMetabolite' && header.match(regexMet) == null)) {
          continue;
        }
        const dbName = dbnameDict[fcomponent]['dbname_map'][header];
        const rawExternalId = utils.cleanExternalId(contentArr[j], dbName);
        if (rawExternalId == '') { //ignore the record whithout any valid externalId
          continue;
        }
        // There might be multiple ids in one externalId item
        externalIdArr = rawExternalId.split(';').map(e => e.trim());

        for (var externalId of externalIdArr) {
          const url_prefix = dbnameDict[fcomponent]['url_map'][header];
          var url = "";
          if ( url_prefix != '' && externalId != '') {
            url = dbnameDict[fcomponent]['url_map'][header] + ':' + externalId;
          }

          const externalDbEntryKey = `${dbName}${externalId}${url}`; // diff url leads to new nodes!
          externalId = dbName === 'ChEBI' ? 'CHEBI:'+externalId : externalId;

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
      }
    }
  } else {
    console.log(`Warning: cannot find external ID file ${filename} in path`, modelDir);
  }

  // write the associated file
  csvWriter = createCsvWriter({
    path: `${outputPath}${fcomponent}ExternalDbs.csv`,
    header: [{ id: `${fcomponent}Id`, title: `${fcomponent}Id` },
              { id: 'externalDbId', title: 'externalDbId' }],
  });
  csvWriter.writeRecords(externalIdDBComponentRel.map(
    (e) => { return { [`${fcomponent}Id`]: e.id, externalDbId: e.externalDbId }; }
  ));
  return extNodeIdTracker;
}

const getUniqueCompartmentlizedMap = (m, hm, uniqueCompartmentalizedMap) => {
  const newID = utils.idfyString2(m.name);
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
}

const getUniqueMetabolite = (m, uniqueCompartmentalizedMap, uniqueMetDict, uniqueMetabolites) => {
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
}

const getModelCypherInstructions = (prefix, dropIndexes, model, version, instructions) => {
// Get cyper instructions for each GEM model
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
  return instructions;
}

const getRemainCyperInstructions = (instructions) => {
// Get the remaining cyper instructions
  `CALL db.index.fulltext.createNodeIndex(
    "fulltext",
    ["CompartmentState", "Compartment", "MetaboliteState", "Metabolite", "CompartmentalizedMetabolite", "SubsystemState", "Subsystem", "ReactionState", "Reaction", "GeneState", "Gene", "PubmedReference"],
    ["id", "name", "letterCode", "alternateName", "synonyms", "description", "formula", "function", "pubMedID", "ec"]);
  `.split('\n').forEach(i => {
    instructions.push(i);
  });
  return instructions;
}
const getGeneIdsFromGeneRule = (geneRule) => {
  let idList = [];
  if (geneRule) {
    idList = geneRule.split(/\s+and\s+|\s+or\s+/).filter(e => e);
  }
  return idList;
}

exports.getGeneIdsFromGeneRule = getGeneIdsFromGeneRule;
exports.createComponentSVGMapFile = createComponentSVGMapFile;
exports.createPMIDFile = createPMIDFile;
exports.extractGeneAnnotation = extractGeneAnnotation;
exports.createComponentExternalDbFile = createComponentExternalDbFile;
exports.getUniqueCompartmentlizedMap = getUniqueCompartmentlizedMap;
exports.getUniqueMetabolite = getUniqueMetabolite;
exports.getModelCypherInstructions = getModelCypherInstructions;
exports.getRemainCyperInstructions = getRemainCyperInstructions;
