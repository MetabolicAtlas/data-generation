const fs = require('fs');
const path = require('path');

const DATA_TYPE_COMPONENTS = {
  transcriptomics: 'gene',
  metabolomics: 'compartmentalizedMetabolite',
};

/*
 * This function transforms data files into files that are ready
 * to be used in the Metabolic Atlas website.
 * Example:
 * `modelDir`: ../data-files/integrated-models/Human-GEM
 * output index file: ./data/dataOverlay/Human-GEM/index.json
 * output data source file: ./data/dataOverlay/Human-GEM/transcriptomics/protein1.mock.tsv
 */
const processDataOverlayFiles = ({ modelDir, outDir, componentIdDict }) => {
  const filesDir = `${modelDir}/dataOverlay`;
  if (!fs.existsSync(filesDir)) {
    return;
  }

  const modelOutDir = getModelOutDir({ modelDir, outDir });

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

  fs.writeFileSync(
    `${modelOutDir}/index.json`,
    JSON.stringify(dataSourcesDict),
    'utf8',
  );

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

      const dataSourceOutDir = `${modelOutDir}/${dt}`;
      if (!fs.existsSync(`${dataSourceOutDir}`)) {
        fs.mkdirSync(`${dataSourceOutDir}`);
      }

      fs.writeFileSync(
        `${dataSourceOutDir}/${filename}`,
        condensedFile,
        'utf8',
      );
    }
  }
};

const getModelOutDir = ({ modelDir, outDir }) => {
  const dataOverlayOutDir = `${outDir}/dataOverlay`;

  if (!fs.existsSync(`${dataOverlayOutDir}`)) {
    fs.mkdirSync(`${dataOverlayOutDir}`);
  }

  const modelFolder = modelDir.split('/').pop();
  const modelOutDir = `${outDir}/dataOverlay/${modelFolder}`;
  if (!fs.existsSync(`${modelOutDir}`)) {
    fs.mkdirSync(`${modelOutDir}`);
  }

  return modelOutDir;
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

  return [header, ...filteredRows].join('\n');
};

module.exports = { processDataOverlayFiles };
