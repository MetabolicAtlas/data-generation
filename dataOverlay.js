const fs = require('fs');
const path = require('path');

const DATA_TYPE_COMPONENTS = {
  transcriptomics: 'gene',
  metabolomics: 'compartmentalizedMetabolite',
};

const processDataOverlayFiles = ({ modelDir, componentIdDict }) => {
  const filesDir = `${modelDir}/dataOverlay`;
  if (!fs.existsSync(filesDir)) {
    return;
  }

  const dataOverlayFiles = {};
  const dataTypes = fs
    .readdirSync(filesDir, { withFileTypes: true })
    .filter((dirent) => dirent.isDirectory())
    .map((dirent) => dirent.name);

  const dataSourcesDict = dataTypes.reduce(
    (obj, dt) => ({
      ...obj,
      [dt]: parseIndexFile(
        fs.readFileSync(`${filesDir}/${dt}/index.tsv`, 'utf8'),
      ),
    }),
    {},
  );

  // TODO: write dataSourcesDict to data folder

  for (const [dt, metadataList] of Object.entries(dataSourcesDict)) {
    const componentType = DATA_TYPE_COMPONENTS[dt];
    const componentIdSet = new Set(Object.keys(componentIdDict[componentType]));

    for (const { filename } of metadataList) {
      const inputFile = fs.readFileSync(
        `${filesDir}/${dt}/${filename}`,
        'utf8',
      );
      const condensedFile = condenseDataSourceFile({
        inputFile,
        componentIdSet,
      });
      // TODO: write condensedFile to data folder
      console.log(condensedFile);
    }
  }
};

const parseIndexFile = (indexFile) => {
  const dataSources = [];
  const [header, ...rows] = indexFile.split('\n').filter(Boolean);
  const keys = header.trim().split('\t').filter(Boolean);

  return rows.map((row) =>
    keys.reduce(
      (obj, key) => ({
        ...obj,
        [key]: row.trim().split('\t').filter(Boolean)[keys.indexOf(key)],
      }),
      {},
    ),
  );
};

const condenseDataSourceFile = ({ inputFile, componentIdSet }) => {
  const [header, ...rows] = inputFile.split('\n').filter(Boolean);

  const filteredRows = rows.filter((row) => {
    const [id] = row.trim().split('\t').filter(Boolean);
    return componentIdSet.has(id);
  });

  return [header, ...filteredRows];
};

module.exports = { processDataOverlayFiles };
