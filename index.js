import fs from 'fs';
import path from 'path';
import * as parser from './parser.js';
import * as utils from './utils.js';
import * as writer from './writer.js';
import * as cypher from './cypher.js';
import { processDataOverlayFiles } from './dataOverlay.js';
import { createTimelineChart } from './timelineGenerator.js';

let extNodeIdTracker = 1;
const humanGeneIdSet = new Set();
const externalIdDBMap = {};
const PMIDSset = new Set();
let instructions = [];
let dropIndexes = false;
let prefix = '';
let outputPath = '';
let outDir = './neo4j';
let inputDir, yamlFile;

const parseModelFiles = (modelDir) => {
  // find the yaml in the folder
  yamlFile = utils.getFile(modelDir, /.*[.](yaml|yml)$/);
  if (!yamlFile) {
    throw new Error('yaml file not found in path ', modelDir);
  }

  const [
    metadata,
    metabolites,
    reactions,
    genes,
    compartments,
    metadataSection,
    model,
    version,
    isHuman,
  ] = parser.getInfoFromYaml(yamlFile);

  prefix = `${model}${version}`;
  outputPath = `${outDir}/${prefix}.`;

  const content = {
    // reformat object as proper key:value objects, rename/add/remove keys
    compartmentalizedMetabolite:
      utils.reformatCompartmentalizedMetaboliteObjets(metabolites.metabolites),
    reaction: utils.reformatReactionObjets(reactions.reactions),
    gene: utils.reformatGeneObjets(genes.genes),
    compartment: utils.reformatCompartmentObjets(compartments.compartments),
  };

  const componentIdDict = utils.getComponentIdDict(content);
  if (isHuman) {
    utils.getHumanGeneIdSet(componentIdDict, humanGeneIdSet);
  }

  processDataOverlayFiles({
    modelDir,
    outDir: './dataOverlay',
    componentIdDict,
  });

  // ========================================================================
  // SVG mapping file
  const svgNodes = [];
  ['compartment', 'subsystem', 'custom'].forEach((component) => {
    const svgRels = parser.getComponentSvgRel(
      content,
      component,
      svgNodes,
      modelDir,
      metadataSection
    );
    writer.writeComponentSvgCSV(svgRels, outputPath, component);
  });

  // write svgMaps file
  writer.writeSvgCSV(svgNodes, outputPath);

  // ========================================================================
  // parse pubmed references
  const [PMIDs, reactionPMID] = parser.getPMIDs(PMIDSset, componentIdDict);
  writer.writePMIDCSV(PMIDs, outputPath); // write pubmedReferences file
  writer.writeReactionPMIDCSV(reactionPMID, outputPath); // write reaction pubmed reference file

  // ========================================================================
  // parse gene annotations
  parser.getGeneAnnotation(componentIdDict, modelDir);

  // ========================================================================
  // parse External IDs files
  const externalIdNodes = [];

  ['reaction', 'metabolite', 'gene'].forEach((component) => {
    let fcomponent = '';
    let externalIdDBComponentRel = [];
    [extNodeIdTracker, fcomponent, externalIdDBComponentRel] =
      parser.getComponentExternalDb(
        externalIdNodes,
        externalIdDBMap,
        extNodeIdTracker,
        component,
        componentIdDict,
        modelDir,
      );
    writer.writeComponentExternalDbCSV(
      externalIdDBComponentRel,
      outputPath,
      fcomponent,
    );
  });

  // write the externalDbs file
  writer.writeExternalDbCSV(externalIdNodes, outputPath);

  // ========================================================================
  // write main nodes relationships files
  writer.writeMetaboliteCompartmentCSV(content, outputPath);

  // ========================================================================
  // write metabolite-compartmentalizedMetabolite relationships
  // generate unique metabolite
  // keep only distinct metabolite (non-compartmentalize) and use the name to generate IDs
  let hm = {};
  const uniqueCompartmentalizedMap = {};
  content.compartmentalizedMetabolite.forEach((m) => {
    utils.getUniqueCompartmentlizedMap(m, hm, uniqueCompartmentalizedMap);
  });

  const uniqueMetDict = {};
  const uniqueMetabolites = [];
  content.compartmentalizedMetabolite.forEach((m) => {
    utils.getUniqueMetabolite(
      m,
      uniqueCompartmentalizedMap,
      uniqueMetDict,
      uniqueMetabolites,
    );
  });

  // create compartmentalizedMetabolite file
  writer.writeMetaboliteCSV(content, outputPath);

  // ========================================================================
  // extract information from metabolite annotation file
  // METABOLITES.tsv has been removed for the format, and this file is actually
  // empty in the old format

  // ========================================================================
  // CM-M relationships
  writer.writeMetaboliteMetaboliteRelCSV(
    content,
    uniqueCompartmentalizedMap,
    outputPath,
  );

  // delete compartmentlizedMetabolites, add unique metabolites
  content.metabolite = uniqueMetabolites;
  delete content.compartmentalizedMetabolite;

  // write reactants-reaction, reaction-products, reaction-genes, reaction-susbsystems relationships files
  const [
    reactionReactantRecords,
    reactionProductRecords,
    reactionGeneRecords,
    reactionSubsystemRecords,
  ] = utils.getReactionRel(content);
  writer.writeRRCSV(reactionReactantRecords, outputPath);
  writer.writeRPCSV(reactionProductRecords, outputPath);
  writer.writeRGCSV(reactionGeneRecords, outputPath);
  writer.writeRSCSV(reactionSubsystemRecords, outputPath);

  // ========================================================================
  // write nodes files
  Object.keys(content).forEach((k) => {
    const elements = content[k];
    writer.writeComponentCSV(content, k, outputPath);
    writer.writeComponentStateCSV(content, k, outputPath);
  });

  // TODO generate instructions more dynamically
  instructions = cypher.getModelCypherInstructions(
    prefix,
    dropIndexes,
    model,
    version,
    instructions,
  );
};

// ========================================================================
// argument parsing
const args = [];
try {
  for (let i = 0; i < process.argv.length; i += 1) {
    if (process.argv[i] === '--reset-db') {
      dropIndexes = true;
    } else {
      args.push(process.argv[i]);
    }
  }
  inputDir = args[2] + '/integrated-models';
} catch {
  console.log('Usage: yarn start input_dir');
  console.log('Usage: yarn start input_dir --reset-db');
}

if (!fs.existsSync(`${outDir}`)) {
  fs.mkdirSync(`${outDir}`);
}

// ========================================================================
// main procedure
try {
  const intputDirFiles = fs.readdirSync(inputDir);
  for (let i = 0; i < intputDirFiles.length; i++) {
    const filePath = path.join(inputDir, intputDirFiles[i]);
    const stat = fs.lstatSync(filePath);
    if (stat.isDirectory() && intputDirFiles[i][0] != '.') {
      parseModelFiles(filePath);
    }
  }
  instructions = cypher.getRemainCypherInstructions(instructions);
  createTimelineChart(inputDir);
} catch (e) {
  if (e.mark) {
    // avoid to print the whole yaml into console
    e.mark.buffer = '';
  }
  console.log(e);
}

// ========================================================================
// write cypher intructions to file
writer.writeCypherFile(instructions, outDir);
