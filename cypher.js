// functions to get cypher instructions
const getModelCypherInstructions = (prefix, dropIndexes, model, version, instructions) => {
// Get cypher instructions for each GEM model
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
  return instructions;
}

const getRemainCypherInstructions = (instructions) => {
// Get the remaining cypher instructions
  `CALL db.index.fulltext.createNodeIndex(
    "fulltext",
    ["CompartmentState", "Compartment", "MetaboliteState", "Metabolite", "CompartmentalizedMetabolite", "SubsystemState", "Subsystem", "ReactionState", "Reaction", "GeneState", "Gene", "PubmedReference"],
    ["id", "name", "letterCode", "alternateName", "synonyms", "description", "formula", "function", "pubMedID", "ec"]);
  `.split('\n').forEach(i => {
    instructions.push(i);
  });
  return instructions;
}

module.exports = {
  getModelCypherInstructions,
  getRemainCypherInstructions,
};
