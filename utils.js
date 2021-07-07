const fs = require('fs'), path = require('path');

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
  if (dbName == 'ChEBI') {
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
}

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
}

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
}

const getReactionRel = (content) => {
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
  return [reactionReactantRecords, reactionProductRecords, reactionGeneRecords, reactionSubsystemRecords];
}

const getGeneIdsFromGeneRule = (geneRule) => {
  let idList = [];
  if (geneRule) {
    idList = geneRule.split(/\s+and\s+|\s+or\s+/)
      .filter(e => e)
      .map(gid => gid.replace(/(\(+|\)+)/, ''));
  }

  // convert to set and back to array
  return [...new Set(idList)];
}

const getUniqueCompartmentlizedMap = (m, hm, uniqueCompartmentalizedMap) => {
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

const getComponentIdDict = (content, isHuman) => {
  const componentIdDict = {}; // store for each type of component the key  Id <-> element
  // use to filter out annotation/external ids for components not in the model and to add missing information
  // extracted from these annotation files such as description, etc...
  Object.keys(content).forEach((k) => {
    componentIdDict[k] = Object.fromEntries(content[k].map(e => [e[`${k}Id`], e]));
  });

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
  return componentIdDict;
}

const getHumanGeneIdSet = (componentIdDict, isHuman, humanGeneIdSet) => {
  // the parameter humanGeneIdSet will be updated
  if (isHuman) {
    Object.keys(componentIdDict.gene).forEach((geneId) => {
      humanGeneIdSet.add(geneId);
    });
  }
}

module.exports = {
  getFile,
  toLabelCase,
  trim,
  cleanExternalId,
  reformatGeneObjets,
  reformatCompartmentObjets,
  reformatCompartmentalizedMetaboliteObjets,
  reformatReactionObjets,
  idfyString,
  idfyString2,
  getReactionRel,
  getGeneIdsFromGeneRule,
  getUniqueCompartmentlizedMap,
  getUniqueMetabolite,
  getComponentIdDict,
  getHumanGeneIdSet,
};

