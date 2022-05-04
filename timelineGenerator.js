import fs from 'fs';
import path from 'path';
import * as d3 from 'd3';
import { JSDOM } from 'jsdom';

const outputPath = 'gemRepository/integratedModelsTimeline.svg';
const d3Curve = d3.line().curve(d3.curveNatural);
const WIDTH = 1500;
let GEMS, HEIGHT, GEM_POSITIONS, BRANCH_POINTS, TIME_SCALE;

const createTimelineChart = async (integratedModelsPath) => {
  GEMS = setGems(integratedModelsPath);
  HEIGHT = getHeight();
  GEM_POSITIONS = getGemPositions();
  BRANCH_POINTS = getBranchPoints();
  TIME_SCALE = createTimeScale();

  const fakeDom = new JSDOM('<!DOCTYPE html><html><body></body></html>');
  const body = (await d3).select(fakeDom.window.document).select('body');

  let svg = body.append('svg').attr('width', WIDTH).attr('height', HEIGHT);

  svg = addBranchPoints(svg);
  svg = addGemVersions(svg);
  svg = addTimeAxis(svg);
  svg = addGemsAxis(svg);

  fs.writeFileSync(outputPath, body.html());
};

const setGems = (integratedModelsPath) => {
  const gems = [];

  const intputDirFiles = fs.readdirSync(integratedModelsPath);
  for (let i = 0; i < intputDirFiles.length; i++) {
    const filePath = path.join(integratedModelsPath, intputDirFiles[i]);
    const stat = fs.lstatSync(filePath);
    if (stat.isDirectory() && intputDirFiles[i][0] != '.') {
      const rawJson = fs.readFileSync(`${filePath}/gemRepository.json`, 'utf8');
      gems.push(JSON.parse(rawJson));
    }
  }

  return gems.map((versions) =>
    versions.filter((v) => v.id.split('.')[2] === '0'),
  );
};

const getHeight = () => GEMS.length * 50 + 100;

const getModelFromId = (id) => id.split('-')[0];

const getGemPositions = () =>
  GEMS.reduce(
    (acc, [v], i) => ({
      ...acc,
      [getModelFromId(v.id)]: HEIGHT - (i + 1) * 50,
    }),
    {},
  );

const getBranchPoints = () => {
  BRANCH_POINTS = {};

  Array.from(GEMS).forEach(([v]) => {
    if (v.externalParentId.length > 0) {
      BRANCH_POINTS[v.id] = v.externalParentId;
    }
  });

  return BRANCH_POINTS;
};

const createTimeScale = () => {
  const allVersionDates = GEMS.flat().map((v) => new Date(v.releaseDate));
  const earliestDate = allVersionDates.sort((a, b) => a - b)[0];
  const timePadding = 60 * 60 * 24 * 30 * 3 * 1000; // 3 months
  const startDate = new Date(earliestDate - timePadding);
  const currentDate = new Date(Date.now() + timePadding);
  return d3.scaleTime().domain([startDate, currentDate]).range([0, WIDTH]);
};

const addTimeAxis = (svg) => {
  const xAxis = d3.axisTop().scale(TIME_SCALE);
  svg
    .append('g')
    .attr('transform', `translate(0, ${HEIGHT})`)
    .style('font-weight', 'bold')
    .call(xAxis);
  return svg;
};

const addGemsAxis = (svg) => {
  const kvFlippedGemPositions = Object.fromEntries(
    Object.entries(GEM_POSITIONS).map(([k, v]) => [v, k]),
  );

  const yScaleLeft = d3.scaleLinear().domain([0, HEIGHT]).range([HEIGHT, 0]);
  const yAxisLeft = d3
    .axisRight(yScaleLeft)
    .tickValues(Object.keys(kvFlippedGemPositions))
    .tickFormat((x) => kvFlippedGemPositions[x]);
  svg
    .append('g')
    .style('font-size', '12px')
    .style('font-weight', 'bold')
    .call(yAxisLeft);

  // const yScaleRight = d3.scaleLinear().domain([0, HEIGHT]).range([HEIGHT, 0]);
  // const yAxisRight = d3
  //   .axisLeft(yScaleRight)
  //   .tickValues(Object.keys(kvFlippedGemPositions))
  //   .tickFormat(x => kvFlippedGemPositions[x]);
  // svg.append('g').attr('transform', `translate(${WIDTH}, 0)`).call(yAxisRight);

  return svg;
};

const addBranchPoints = (svg) => {
  const allGemVersions = GEMS.flat();

  Object.entries(BRANCH_POINTS).forEach(([childId, parentIds]) => {
    const childVersion = allGemVersions.find((v) => v.id === childId);
    const childModel = getModelFromId(childVersion.id);

    parentIds.forEach(({ id, citLink }) => {
      const parentVersion = allGemVersions.find((v) => v.id === id);
      const parentModel = getModelFromId(parentVersion.id);

      const curve = [
        [
          TIME_SCALE(new Date(parentVersion.releaseDate)),
          HEIGHT - GEM_POSITIONS[parentModel],
        ],
        [
          TIME_SCALE(new Date(parentVersion.releaseDate)) +
            (TIME_SCALE(new Date(childVersion.releaseDate)) -
              TIME_SCALE(new Date(parentVersion.releaseDate))) /
              5,
          HEIGHT -
            GEM_POSITIONS[childModel] -
            (HEIGHT -
              GEM_POSITIONS[childModel] -
              (HEIGHT - GEM_POSITIONS[parentModel])) /
              5,
        ],
        [
          TIME_SCALE(new Date(childVersion.releaseDate)),
          HEIGHT - GEM_POSITIONS[childModel],
        ],
      ];

      svg
        .append('path')
        .attr('d', d3Curve(curve))
        .attr('stroke', 'lightgray')
        .attr('fill', 'none');
    });
  });

  return svg;
};

const addGemVersions = (svg) => {
  GEMS.forEach((versions) => {
    const model = getModelFromId(versions[0].id);
    const y = HEIGHT - GEM_POSITIONS[model];

    // add connecting lines
    svg
      .append('g')
      .selectAll('line')
      .data(versions.slice(0, -1).map((v) => new Date(v.releaseDate)))
      .enter()
      .append('line')
      .style('stroke', 'lightgray')
      .attr('x1', (d) => TIME_SCALE(d))
      .attr('y1', y)
      .attr('x2', (_d, i) => {
        const nextD = new Date(versions[i + 1].releaseDate);
        return TIME_SCALE(nextD);
      })
      .attr('y2', y);

    //  add circles
    svg
      .append('g')
      .selectAll('circle')
      .data(versions.map((v) => new Date(v.releaseDate)))
      .join('circle')
      .classed('circle', true)
      .attr('r', (_d, i) => {
        const idParts = versions[i].id.split('.');
        const isMajor = idParts[1] === '0' && idParts[2] === '0';
        return isMajor ? 24 : 18;
      })
      .attr('cy', y)
      .attr('cx', (d) => TIME_SCALE(d))
      .attr('data-model', (_d, i) => versions[i].id.match(/\w+-\w+/)[0])
      .attr('data-version', (_d, i) => versions[i].id.split('-')[2])
      .attr('data-release-date', (_d, i) => versions[i].releaseDate)
      .attr('data-release-link', (_d, i) => versions[i].releaseLink)
      .attr('data-pmid', (_d, i) => versions[i].PMID)
      .attr('data-external-parent-ids', (_d, i) =>
        JSON.stringify(versions[i].externalParentId),
      );

    svg
      .append('g')
      .selectAll('text')
      .data(versions.map((v) => new Date(v.releaseDate)))
      .join('text')
      .classed('label', true)
      .style('fill', 'white')
      .style('font-size', '10px')
      .style('font-weight', 'bold')
      .attr('y', y + 3)
      .attr(
        'x',
        (d, i) => TIME_SCALE(d) - (versions[i].id.split('.')[1] > 9 ? 5 : 2),
      )
      .attr('transform', 'translate(-10, 0)')
      .text((_d, i) => versions[i].id.split('-')[2]);
  });

  return svg;
};

export { createTimelineChart };
