const fs = require('fs'), path = require('path');
const yaml = require('js-yaml');
const func = require('./func.js');

const createCsvWriter = require("csv-writer").createObjectCsvWriter;
let csvWriter = null;

const idfyString = s => s.toLowerCase().replace(/[^a-z0-9]/g, '_').replace(/_+/g, '_').replace(/^_|_$/, ''); // for subsystems, compartments etc..
const idfyString2 = s => s.toLowerCase().replace(/[^a-z0-9]/g, '_').replace(/_+/g, '_'); // to generate compartmentalizedMetabolite ID from their name

const mergedObjects = data => data.reduce((acc, item) => {
  const [key, value] = Object.entries(item)[0];
    return {
    ...acc,
        [key]: value,
    };
}, {});


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
      reversible: r.lower_bound === -1000,
      ec: r.eccodes,
      references: r.references,
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
    const filename = `${component}SVG.tsv.plain`;
    const mappingFile = func.getFile(modelDir, filename);
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
            [`${component}Id`]: idfyString(componentName),
            svgMapId: mapFilename.split('.')[0],
          });
        }
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
    csvWriter.writeRecords(svgRels);
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

  // extract information from gene annotation file
  const geneAnnoFile = func.getFile(modelDir, /genes-new[.]tsv$/);
  if (!geneAnnoFile) {
    console.log("Warning: cannot find gene annotation file genes.tsv in path", modelDir);
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
      const [ geneId, geneENSTID, geneENSPID, geneUniProtID, name, geneEntrezID, alternateName, synonyms] = lines[i].split('\t').map(e => func.trim(e, '"'));
      if (geneId in componentIdDict.gene) { //only keep the ones in the model
        const gene = componentIdDict.gene[geneId];
        Object.assign(gene, { name, alternateName, synonyms, function: thefunction }); // other props are not in the db design, TODO remove them?
      }
    }
  }

  // extract description subsystem annotation file
  // TODO or remove annotation file

  // ========================================================================
  // parse External IDs files

  const { dbnameDict } = require('./var');
  const externalIdNodes = [];

  ['reaction', 'metabolite', 'gene', 'subsystem'].forEach((component) => {
    const externalIdDBComponentRel = [];
    const filename = `${component}s-new.tsv`;
    const extIDFile = func.getFile(modelDir, filename);
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
          contentArr = lines[i].split('\t').map(e => func.trim(e, '"'));
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
          const rawExternalId = func.cleanExternalId(contentArr[j], dbName);
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
      reactionSubsystemRecords.push({ reactionId: r.reactionId, subsystemId: idfyString(name) });
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

