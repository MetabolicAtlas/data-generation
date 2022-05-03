import fs from 'fs';
import path from 'path';
import { createObjectCsvWriter } from 'csv-writer';
import * as utils from './utils.js';

const createCsvWriter = createObjectCsvWriter;

const writeComponentSvgCSV = (svgRels, outputPath, component) => {
  // Write CSV files with componentId and SvgMapId
  const csvWriter = createCsvWriter({
    path: `${outputPath}${component}SvgMaps.csv`,
    header: [
      { id: `${component}Id`, title: `${component}Id` },
      { id: 'svgMapId', title: 'svgMapId' },
    ],
  });
  csvWriter.writeRecords(svgRels);
};

const writePMIDCSV = (PMIDs, outputPath) => {
  // write the file containing a list of pubmed IDs
  const csvWriter = createCsvWriter({
    path: `${outputPath}pubmedReferences.csv`,
    header: [{ id: 'id', title: 'id' }],
  });
  csvWriter.writeRecords(
    PMIDs.map((id) => {
      return { id };
    }),
  );
};
const writeReactionPMIDCSV = (reactionPMID, outputPath) => {
  // write the file containing reaction IDs and its corresponding pubmed IDs
  const csvWriter = createCsvWriter({
    path: `${outputPath}reactionPubmedReferences.csv`,
    header: [
      { id: 'reactionId', title: 'reactionId' },
      { id: 'pubmedReferenceId', title: 'pubmedReferenceId' },
    ],
  });
  csvWriter.writeRecords(reactionPMID);
};

const writeComponentExternalDbCSV = (
  externalIdDBComponentRel,
  outputPath,
  fcomponent,
) => {
  // Write the file containing externalDb IDs for each component
  const csvWriter = createCsvWriter({
    path: `${outputPath}${fcomponent}ExternalDbs.csv`,
    header: [
      { id: `${fcomponent}Id`, title: `${fcomponent}Id` },
      { id: 'externalDbId', title: 'externalDbId' },
    ],
  });
  csvWriter.writeRecords(
    externalIdDBComponentRel.map((e) => {
      return { [`${fcomponent}Id`]: e.id, externalDbId: e.externalDbId };
    }),
  );
};

const writeSvgCSV = (svgNodes, outputPath) => {
  // Write the file containing a list of SvgMap files
  const csvWriter = createCsvWriter({
    path: `${outputPath}svgMaps.csv`,
    header: svgNodes.length
      ? Object.keys(svgNodes[0]).map((k) => Object({ id: k, title: k }))
      : '',
  });
  csvWriter.writeRecords(svgNodes);
};

const writeExternalDbCSV = (externalIdNodes, outputPath) => {
  // Write the file containing all externalDb info
  const csvWriter = createCsvWriter({
    path: `${outputPath}externalDbs.csv`,
    header: [
      { id: 'id', title: 'id' },
      { id: 'dbName', title: 'dbName' },
      { id: 'externalId', title: 'externalId' },
      { id: 'url', title: 'url' },
    ],
  });
  csvWriter.writeRecords(externalIdNodes);
};

const writeMetaboliteCompartmentCSV = (content, outputPath) => {
  const compartmentLetterToIdMap = content.compartment.reduce((entries, c) => {
    // return a map which can get the compartment ID from the compartment letter
    return {
      ...entries,
      [c.letterCode]: c.compartmentId,
    };
  }, {});

  const csvWriter = createCsvWriter({
    path: `${outputPath}compartmentalizedMetaboliteCompartments.csv`,
    header: [
      {
        id: 'compartmentalizedMetaboliteId',
        title: 'compartmentalizedMetaboliteId',
      },
      { id: 'compartmentId', title: 'compartmentId' },
    ],
  });

  csvWriter.writeRecords(
    content.compartmentalizedMetabolite.map((e) => {
      return {
        compartmentalizedMetaboliteId: e.compartmentalizedMetaboliteId,
        compartmentId: compartmentLetterToIdMap[e.compartment],
      };
    }),
  );
};

const writeMetaboliteCSV = (content, outputPath) => {
  const csvWriter = createCsvWriter({
    path: `${outputPath}compartmentalizedMetabolites.csv`,
    header: [{ id: 'id', title: 'id' }],
  });

  csvWriter.writeRecords(
    content.compartmentalizedMetabolite.map((e) => {
      return { id: e.compartmentalizedMetaboliteId };
    }),
  );
};

const writeMetaboliteMetaboliteRelCSV = (
  content,
  uniqueCompartmentalizedMap,
  outputPath,
) => {
  const csvWriter = createCsvWriter({
    path: `${outputPath}compartmentalizedMetaboliteMetabolites.csv`,
    header: [
      {
        id: 'compartmentalizedMetaboliteId',
        title: 'compartmentalizedMetaboliteId',
      },
      { id: 'metaboliteId', title: 'metaboliteId' },
    ],
  });

  csvWriter.writeRecords(
    content.compartmentalizedMetabolite.map((e) => {
      return {
        compartmentalizedMetaboliteId: e.compartmentalizedMetaboliteId,
        metaboliteId:
          uniqueCompartmentalizedMap[e.compartmentalizedMetaboliteId],
      };
    }),
  );
};

const writeRRCSV = (reactionReactantRecords, outputPath) => {
  // write reaction-reactants relationship files
  const csvWriterRR = createCsvWriter({
    path: `${outputPath}compartmentalizedMetaboliteReactions.csv`,
    header: [
      {
        id: 'compartmentalizedMetaboliteId',
        title: 'compartmentalizedMetaboliteId',
      },
      { id: 'reactionId', title: 'reactionId' },
      { id: 'stoichiometry', title: 'stoichiometry' },
    ],
  });
  csvWriterRR.writeRecords(reactionReactantRecords);
};
const writeRPCSV = (reactionProductRecords, outputPath) => {
  // write reaction-products relationship files
  const csvWriterRP = createCsvWriter({
    path: `${outputPath}reactionCompartmentalizedMetabolites.csv`,
    header: [
      { id: 'reactionId', title: 'reactionId' },
      {
        id: 'compartmentalizedMetaboliteId',
        title: 'compartmentalizedMetaboliteId',
      },
      { id: 'stoichiometry', title: 'stoichiometry' },
    ],
  });
  csvWriterRP.writeRecords(reactionProductRecords);
};
const writeRGCSV = (reactionGeneRecords, outputPath) => {
  // write reaction-genes relationship files
  const csvWriterRG = createCsvWriter({
    path: `${outputPath}reactionGenes.csv`,
    header: [
      { id: 'reactionId', title: 'reactionId' },
      { id: 'geneId', title: 'geneId' },
    ],
  });
  csvWriterRG.writeRecords(reactionGeneRecords);
};
const writeRSCSV = (reactionSubsystemRecords, outputPath) => {
  // write reaction-subsystems relationship files
  const csvWriterRS = createCsvWriter({
    path: `${outputPath}reactionSubsystems.csv`,
    header: [
      { id: 'reactionId', title: 'reactionId' },
      { id: 'subsystemId', title: 'subsystemId' },
    ],
  });
  csvWriterRS.writeRecords(reactionSubsystemRecords);
};

const writeComponentCSV = (content, k, outputPath) => {
  const elements = content[k];
  const csvWriter = createCsvWriter({
    path: `${outputPath}${k}s.csv`,
    header: [Object({ id: 'id', title: 'id' })],
  });
  csvWriter.writeRecords(elements.map((e) => Object({ id: e[`${k}Id`] })));
};

const writeComponentStateCSV = (content, k, outputPath) => {
  const elements = content[k];
  const csvWriter = createCsvWriter({
    path: `${outputPath}${k}States.csv`,
    header: Object.keys(elements[0])
      // ignore some keys 'metabolites', 'subsystems' are in reactions, 'compartment' is in metabolite
      .filter((k) => !['metabolites', 'subsystems', 'compartment'].includes(k))
      .map((k) => Object({ id: k, title: k })),
  });
  // destructure object to remove the keys
  csvWriter.writeRecords(
    elements.map(({ subsystems, metabolites, compartment, ...e }) => e),
  );
};

const writeCypherFile = (instructions, outDir) => {
  fs.writeFileSync(`${outDir}/import.cypher`, instructions.join('\n'), 'utf8');
};

export {
  writeComponentSvgCSV,
  writePMIDCSV,
  writeReactionPMIDCSV,
  writeComponentExternalDbCSV,
  writeSvgCSV,
  writeExternalDbCSV,
  writeMetaboliteCompartmentCSV,
  writeMetaboliteCSV,
  writeMetaboliteMetaboliteRelCSV,
  writeRRCSV,
  writeRPCSV,
  writeRGCSV,
  writeRSCSV,
  writeComponentCSV,
  writeComponentStateCSV,
  writeCypherFile,
};
