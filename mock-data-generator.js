const createCsvWriter = require("csv-writer").createObjectCsvWriter;

// ENUM
const fieldTypes = {
  STRING: "string",
  INT: "int",
  FLOAT: "float",
  BOOL: "bool",
  STRING_LIST: "string_list"
};

// UTIL
const getRandInt = max => Math.floor(Math.random() * (max || 1000));

const getRandFloat = () => (Math.random() * (10 - 0) + 0.0001).toFixed(3);

const getRandString = () =>
  Math.random()
    .toString(36)
    .substr(2, 5);

const getRandBool = () => Math.random() >= 0.5;

const shuffle = a => {
  for (let i = a.length - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1));
    [a[i], a[j]] = [a[j], a[i]];
  }
  return a;
};

const generateFieldValue = field => {
  const fieldName = field[0];
  const fieldType = field[1];
  let fieldValue;

  switch (fieldType) {
    case fieldTypes.STRING:
      fieldValue = getRandString();
      break;
    case fieldTypes.INT:
      fieldValue = getRandInt();
      break;
    case fieldTypes.FLOAT:
      fieldValue = getRandFloat();
      break;
    case fieldTypes.BOOL:
      fieldValue = getRandBool();
      break;
    case fieldTypes.STRING_LIST:
      fieldValue = [getRandString(), getRandString(), getRandString()];
      break;
  }

  return fieldValue;
};

const makeHeader = fields =>
  fields.map(field => ({ id: field[0], title: field[0] }));

const makeRecords = (fields, total) =>
  [...Array(total).keys()].map(() => {
    const record = {};

    for (let field of fields) {
      const fieldName = field[0];
      record[fieldName] = generateFieldValue(field);
    }

    return record;
  });

const getRecordIdForRel = (records, isUnique) =>
  isUnique ? records.shift() : records[getRandInt(records.length)];

const makeRels = (fields, recordsA, recordsB, isUnique, relFields) => {
  if (isUnique) {
    shuffle(recordsB);
  }

  return recordsA.map(rA => ({
    [fields[0][0]]: rA.id,
    [fields[1][0]]: getRecordIdForRel(recordsB, isUnique).id,
    ...relFields.reduce(
      (obj, rf) => ({ ...obj, [rf[0]]: generateFieldValue(rf) }),
      {}
    )
  }));
};

let csvWriter;

const createCsvFile = ({ header, records, filename }) => {
  csvWriter = createCsvWriter({
    path: `./data/${filename}.csv`,
    header
  });

  csvWriter.writeRecords(records).then(() => {
    console.log(`${filename} generated.`);
  });
};

const buildFields = fields =>
  fields.reduce((res, field) => {
    const fieldName = field[0];
    const fieldType = field[1];

    let val = `csvLine.${fieldName}`;
    if (fieldType === fieldTypes.INT) {
      val = `toInteger(${val})`;
    } else if (fieldType === fieldTypes.BOOL) {
      val = `toBoolean(${val})`;
    } else if (fieldType === fieldTypes.FLOAT) {
      val = `toFloat(${val})`;
    }

    return { ...res, [fieldName]: val };
  }, {});

const buildCypherNodesWithStatesInstructions = ({ nodeType, fields }) => {
  const hasModelLabel = SCHEMA.nodeTypes[nodeType].hasModelLabel;
  const stateFields = SCHEMA.nodeTypes[nodeType].stateFields;

  const nodeLabel = nodeType.charAt(0).toUpperCase() + nodeType.slice(1);
  const nodeFilename = `${nodeType}s`;
  const relId = `${nodeType}Id`;

  const nodeFields = buildFields(fields);

  let instructions = `
CREATE INDEX FOR (n:${nodeLabel}) ON (n.id);
LOAD CSV WITH HEADERS FROM "file:///${nodeFilename}.csv" AS csvLine
CREATE (n:${nodeLabel}${hasModelLabel ? ":" + MODEL : ""} ${JSON.stringify(
    nodeFields
  ).replace(/['"]+/g, "")});`;

  if (!!stateFields) {
    const stateNodeLabel = `${nodeLabel}State`;
    const stateNodeFilename = `${nodeType}States`;
    const stateNodeFields = buildFields(stateFields.slice(1));

    instructions += `
LOAD CSV WITH HEADERS FROM "file:///${stateNodeFilename}.csv" AS csvLine
MATCH (n:${nodeLabel} {id: csvLine.${relId}})
CREATE (ns:${stateNodeLabel} ${JSON.stringify(stateNodeFields).replace(
      /['"]+/g,
      ""
    )})
CREATE (n)-[:V1]->(ns);`;
  }

  instructions += "\n";

  return instructions;
};

const buildCypherRelsInstructions = ({ nodeType1, nodeType2, fields }) => {
  const nodeLabel1 = nodeType1.charAt(0).toUpperCase() + nodeType1.slice(1);
  const nodeLabel2 = nodeType2.charAt(0).toUpperCase() + nodeType2.slice(1);
  const relId1 = `${nodeType1}Id`;
  const relId2 = `${nodeType2}Id`;
  const fileName = `${nodeType1}${nodeLabel2}s`;

  let fieldsInstructions = "";
  if (fields) {
    const relFields = buildFields(fields);
    fieldsInstructions = ` ${JSON.stringify(relFields).replace(/['"]+/g, "")}`;
  }

  return `
LOAD CSV WITH HEADERS FROM "file:///${fileName}.csv" AS csvLine
MATCH (n1:${nodeLabel1} {id: csvLine.${relId1}}),(n2:${nodeLabel2} {id: csvLine.${relId2}})
CREATE (n1)-[:V1${fieldsInstructions}]->(n2);
`;
};

// SCHEMA

const SCHEMA = {
  nodeTypes: {
    metabolite: {
      hasModelLabel: true,
      quantity: 2000,
      fields: [["id", fieldTypes.STRING]],
      stateFields: [
        ["metaboliteId", fieldTypes.STRING],
        ["name", fieldTypes.STRING],
        ["alternateName", fieldTypes.STRING],
        ["synonyms", fieldTypes.STRING_LIST],
        ["description", fieldTypes.STRING],
        ["formula", fieldTypes.STRING],
        ["charge", fieldTypes.INT],
        ["isCurrency", fieldTypes.BOOL]
      ]
    },
    compartment: {
      hasModelLabel: true,
      quantity: 10,
      fields: [["id", fieldTypes.STRING]],
      stateFields: [
        ["compartmentId", fieldTypes.STRING],
        ["name", fieldTypes.STRING],
        ["letterCode", fieldTypes.STRING]
      ]
    },
    compartmentalizedMetabolite: {
      hasModelLabel: true,
      quantity: 6000,
      fields: [["id", fieldTypes.STRING]]
    },
    reaction: {
      quantity: 15000,
      hasModelLabel: true,
      fields: [["id", fieldTypes.STRING]],
      stateFields: [
        ["reactionId", fieldTypes.STRING],
        ["name", fieldTypes.STRING],
        ["reversible", fieldTypes.BOOL],
        ["lowerBound", fieldTypes.INT],
        ["upperBound", fieldTypes.INT],
        ["geneRule", fieldTypes.STRING],
        ["ec", fieldTypes.STRING]
      ]
    },
    gene: {
      quantity: 4000,
      hasModelLabel: true,
      fields: [["id", fieldTypes.STRING]],
      stateFields: [
        ["geneId", fieldTypes.STRING],
        ["name", fieldTypes.STRING],
        ["alternateName", fieldTypes.STRING],
        ["synonyms", fieldTypes.STRING_LIST],
        ["function", fieldTypes.STRING]
      ]
    },
    pubmedReference: {
      quantity: 15000,
      fields: [["id", fieldTypes.STRING], ["pubmedId", fieldTypes.STRING]]
    },
    subsystem: {
      quantity: 150,
      hasModelLabel: true,
      quantity: 50,
      fields: [["id", fieldTypes.STRING]],
      stateFields: [
        ["subsystemId", fieldTypes.STRING],
        ["name", fieldTypes.STRING]
      ]
    },
    svgMap: {
      quantity: 60,
      fields: [
        ["id", fieldTypes.STRING],
        ["filename", fieldTypes.STRING],
        ["customName", fieldTypes.STRING]
      ]
    },
    externalDb: {
      quantity: 50000,
      fields: [
        ["id", fieldTypes.STRING],
        ["dbName", fieldTypes.STRING],
        ["externalId", fieldTypes.STRING],
        ["url", fieldTypes.STRING]
      ]
    }
  },
  relationships: [
    {
      node1: "compartmentalizedMetabolite",
      node2: "metabolite"
    },
    {
      node1: "compartmentalizedMetabolite",
      node2: "compartment"
    },
    {
      node1: "compartmentalizedMetabolite",
      node2: "reaction",
      relFields: [["stoichiometry", fieldTypes.FLOAT]]
    },
    {
      node1: "reaction",
      node2: "compartmentalizedMetabolite",
      relFields: [["stoichiometry", fieldTypes.FLOAT]]
    },
    {
      node1: "reaction",
      node2: "gene"
    },
    {
      node1: "reaction",
      node2: "pubmedReference"
    },
    {
      node1: "reaction",
      node2: "subsystem"
    },
    {
      node1: "compartment",
      node2: "svgMap",
      isUnique: true
    },
    {
      node1: "subsystem",
      node2: "svgMap",
      isUnique: true
    },
    {
      node1: "metabolite",
      node2: "externalDb"
    },
    {
      node1: "subsystem",
      node2: "externalDb"
    },
    {
      node1: "reaction",
      node2: "externalDb"
    },
    {
      node1: "gene",
      node2: "externalDb"
    }
  ]
};

// DATA GENERATION

const generateData = () => {
  const data = {};
  const nodeTypes = Object.keys(SCHEMA.nodeTypes);
  for (let nodeType of nodeTypes) {
    const fields = SCHEMA.nodeTypes[nodeType].fields;
    const stateFields = SCHEMA.nodeTypes[nodeType].stateFields;
    const quantity = SCHEMA.nodeTypes[nodeType].quantity || 100;
    let records;

    if (stateFields) {
      const stateRecords = makeRecords(stateFields, quantity);
      createCsvFile({
        header: makeHeader(stateFields),
        records: stateRecords,
        filename: `${nodeType}States`
      });

      records = stateRecords.map(r => ({
        id: r[stateFields[0][0]]
      }));
    } else {
      records = makeRecords(fields, quantity);
    }

    createCsvFile({
      header: makeHeader(fields),
      records,
      filename: `${nodeType}s`
    });

    data[nodeType] = records;
  }

  for (let rel of SCHEMA.relationships) {
    const { node1, node2, isUnique, relFields = [] } = rel;

    const fields = [
      [`${node1}Id`, fieldTypes.STRING],
      [`${node2}Id`, fieldTypes.STRING]
    ];

    const header = makeHeader([...fields, ...relFields]);

    const records = makeRels(
      fields,
      data[node1],
      data[node2],
      isUnique,
      relFields
    );

    const filename = `${node1}${node2.charAt(0).toUpperCase()}${node2.slice(
      1
    )}s`;

    createCsvFile({
      header,
      records,
      filename
    });
  }

  const buildCypherInstructions = () => {
    let instructions = "";

    for (let nodeType of nodeTypes) {
      instructions += buildCypherNodesWithStatesInstructions({
        nodeType,
        fields: SCHEMA.nodeTypes[nodeType].fields
      });
    }

    for (let rel of SCHEMA.relationships) {
      const { node1, node2, relFields } = rel;

      instructions += buildCypherRelsInstructions({
        nodeType1: node1,
        nodeType2: node2,
        fields: relFields
      });
    }

    return instructions;
  };

  console.log(
    "\n*** CYPHER INSTRUCTIONS START ***\n",
    buildCypherInstructions(),
    "\n*** CYPHER INSTRUCTIONS END ***\n"
  );
};

let match;
const regex = /model=(.+)/;

try {
  match = process.argv.slice(2)[0].match(regex);
  if (!match) {
    throw new Error();
  }
} catch {
  console.log("❗️  Please provide a model, e.g. 'yarn start model=Human'");
  return;
}

const MODEL = match[1];
generateData();
