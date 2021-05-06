const fs = require('fs'), path = require('path');
const createCsvWriter = require("csv-writer").createObjectCsvWriter;
const { dbnameDict } = require('./var');

const trim = (x, characters=" \t\w") => {
  var start = 0;
  while (characters.indexOf(x[start]) >= 0) {
    start += 1;
  }
  var end = x.length - 1;
  while (characters.indexOf(x[end]) >= 0) {
    end -= 1;
  }
  return x.substr(start, end - start + 1);
}

const cleanExternalId = (rawExternalId, dbName) => {
  // clean rawExternalId
  rawExternalId = trim(rawExternalId.trim(), '"');
  if (dbName == 'MA'){
    rawExternalId = rawExternalId.replace(/^MA-/, '');
  } else if (dbName == 'ChEBI') {
    rawExternalId = rawExternalId.replace(/^CHEBI:/, '');
  } else if (dbName == 'Rhea' || dbName == 'RheaMaster') {
    rawExternalId = rawExternalId.replace(/^RHEA:/, '');
  }
  return rawExternalId;
}

const getFile = (dirPath, regexpOrString) => {
  if (!fs.existsSync(dirPath)){
    console.log("no dir ", dirPath);
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
}

const getGeneIdsFromGeneRule = (geneRule) => {
  let idList = [];
  if (geneRule) {
    idList = geneRule.split(/\s+and\s+|\s+or\s+/).filter(e => e);
  }
  return idList;
}

const toLabelCase = (modelName) => {
  return modelName.replace('-', ' ').split(/\s/g).map(word => `${word[0].toUpperCase()}${word.slice(1).toLowerCase()}`).join('');
}

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

const createComponentSVGMapFile = (component, outputPath, svgNodes, modelDir) => {
  const filename = `${component}SVG.tsv.plain`;
  const mappingFile = getFile(modelDir, filename);
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
  const geneAnnoFile = getFile(modelDir, /genes-new[.]tsv$/);
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
      const [ geneId, geneENSTID, geneENSPID, geneUniProtID, name, geneEntrezID, alternateName, synonyms] = lines[i].split('\t').map(e => trim(e, '"'));
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
  const extIDFile = getFile(modelDir, filename);
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
        contentArr = lines[i].split('\t').map(e => trim(e, '"'));
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
        const rawExternalId = cleanExternalId(contentArr[j], dbName);
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

exports.getFile = getFile;
exports.toLabelCase = toLabelCase;
exports.getGeneIdsFromGeneRule = getGeneIdsFromGeneRule;
exports.trim = trim;
exports.cleanExternalId = cleanExternalId;
exports.reformatGeneObjets = reformatGeneObjets;
exports.reformatCompartmentObjets = reformatCompartmentObjets;
exports.reformatCompartmentalizedMetaboliteObjets = reformatCompartmentalizedMetaboliteObjets;
exports.reformatReactionObjets = reformatReactionObjets;
exports.idfyString = idfyString;
exports.idfyString2 = idfyString2;
exports.createComponentSVGMapFile = createComponentSVGMapFile;
exports.createPMIDFile = createPMIDFile;
exports.extractGeneAnnotation = extractGeneAnnotation;
exports.createComponentExternalDbFile = createComponentExternalDbFile;

