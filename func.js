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
}
