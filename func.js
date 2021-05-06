const fs = require('fs'), path = require('path');
module.exports = {
  trim: function (x, characters=" \t\w") {
    var start = 0;
    while (characters.indexOf(x[start]) >= 0) {
      start += 1;
    }
    var end = x.length - 1;
    while (characters.indexOf(x[end]) >= 0) {
      end -= 1;
    }
    return x.substr(start, end - start + 1);
  },

  cleanExternalId: function(rawExternalId, dbName) {
    // clean rawExternalId
    rawExternalId = this.trim(rawExternalId.trim(), '"');
    if (dbName == 'MA'){
      rawExternalId = rawExternalId.replace(/^MA-/, '');
    } else if (dbName == 'ChEBI') {
      rawExternalId = rawExternalId.replace(/^CHEBI:/, '');
    } else if (dbName == 'Rhea' || dbName == 'RheaMaster') {
      rawExternalId = rawExternalId.replace(/^RHEA:/, '');
    }
    return rawExternalId;
  },

  getFile: function(dirPath, regexpOrString) {
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
  },

  getGeneIdsFromGeneRule: function (geneRule) {
    let idList = [];
    if (geneRule) {
      idList = geneRule.split(/\s+and\s+|\s+or\s+/).filter(e => e);
    }
    return idList;
  },

  toLabelCase: function (modelName) {
    return modelName.replace('-', ' ').split(/\s/g).map(word => `${word[0].toUpperCase()}${word.slice(1).toLowerCase()}`).join('');
  },

};
