function adjacencyGraph() {
  var genWidth = 800;
  var genHeight = 500;

  /* ------ Visualization set up ------ */
  var adjacency = d3.select("div#adj")
    .insert("svg:svg", ":first-child")
    .attr("width", genWidth)
    .attr("height", genHeight);

  var edgeColor = d3.scaleLinear().domain([1, 20])
    .range(["whitesmoke", "silver"]);

  var simulation = d3.forceSimulation()
    .force("link", d3.forceLink()
      .id(function (d) { return d.id; })
      .strength(function (d) {
        if (d.className == 'self') return 0.09;
        return 0.02;
      }))
    .force('limit', d3.forceLimit()
      .x0(10)
      .x1(genWidth - 10)
      .y0(10)
      .y1(genHeight - 10))
    .force("charge", d3.forceManyBody())
    .force("center", d3.forceCenter(genWidth / 2, genHeight / 2));

  adjacency.append("g")
    .attr("class", "links");

  adjacency.append("g")
    .attr("class", "nodes");

  var dragNodes = d3.drag();

  function dragstarted(d) {
    if (!d3.event.active) simulation.alphaTarget(0.3).restart();
    simulation.restart();
  }

  function dragged(d) {
    d.fx = d3.event.x;
    d.fy = d3.event.y;
  }

  function dragended(d) {
    if (!d3.event.active) simulation.alphaTarget(0);
    d.fx = null;
    d.fy = null;
    simulation.alphaTarget(0.1);
  }

  function remove() {
    adjacency.selectAll("g").selectAll('*').remove();
  }

  /* ------ Draw simulation ------ */
  function draw(data) {
    /* Create links */
    var linksGroup = adjacency.select("g.links").selectAll("line")
      .data(data.links, function (d) {
        return `${d.source.id}-${d.target.id}`;
      });

    var linksEnter = linksGroup
      .enter();

    linksEnter
      .append("line")
      .attr("stroke-width", function (d) { return 2 * Math.sqrt(d.count); })
      .attr("stroke", function (d) {
        if (d.className == 'self') {
          return d3.color(d.color).brighter();
        }
        return edgeColor(d.count);
      })

    var proteins = linksEnter
      .filter(function (d) {
        return d.className == 'self';
      })
      .append("svg:text")
      .attr("class", "textClass")
      .text(function (d) { return d.protName; });

    var links = adjacency.select("g.links").selectAll('line');
    var proteins = adjacency.select("g.links").selectAll('text');

    linksGroup.exit().remove();

    /* Create nodes */
    var nodeGroup = adjacency.select("g.nodes").selectAll(".node")
      .data(data.nodes, function (d) {
        return d.id;
      });

    var nodeEnter = nodeGroup
      .enter()
      .append("g")
      .attr("class", "node")
      .call(d3.drag()
        .on("start", dragstarted)
        .on("drag", dragged)
        .on("end", dragended))

    nodeEnter
      .append(function (d) {
        const r = 1.2 * (Math.sqrt(d.size));
        if (d.id.indexOf("start") >= 0) {
          const circle = document.createElementNS("http://www.w3.org/2000/svg", "circle");
          circle.setAttribute("r", r);
          return circle;
        } else {
          const hex = document.createElementNS("http://www.w3.org/2000/svg", "path");
          const a = r / 2 * Math.sqrt(3);
          hex.setAttribute("d", `M 0,${-r} L ${a},${-r / 2} L ${a},${r / 2} L 0,${r} L ${-a},${r / 2} L ${-a},${-r / 2} L 0,${-r} Z`);
          hex.setAttribute("stroke", "#fff");
          hex.setAttribute("stroke-width", "1.5")
          return hex;
        }
      })
      .attr("fill", function (d) { return d.color; });

    nodeEnter
      .append("svg:text")
    // .attr("class", "textClass")
    // .text(function (d) { return d.id; });

    nodeGroup = nodeEnter.merge(nodeGroup);

    nodeGroup.exit().remove();

    /* Create simulation */
    simulation
      .nodes(data.nodes)
      .on("tick", function () {
        links
          .attr("x1", function (d) { return d.source.x; })
          .attr("y1", function (d) { return d.source.y; })
          .attr("x2", function (d) { return d.target.x; })
          .attr("y2", function (d) { return d.target.y; });

        proteins
          .attr("x", function (d) {
            return ((d.source.x + d.target.x) / 2);
          })
          .attr("y", function (d) {
            return ((d.source.y + d.target.y)) / 2;
          });

        nodeGroup
          .attr("transform", function (d) {
            return "translate(" + d.x + "," + d.y + ")";
          })
      });

    simulation.force("link")
      .links(data.links);
  }

  /* ------ Slider ------ */
  function createSlider(data) {
    var x = d3.scaleLinear()
      .domain([1, 10])
      .range([1, genWidth - 100])
      .clamp(true);

    d3.select("#slider").select("svg").remove();

    var svg2 = d3.select("div#slider")
      .append("svg:svg")
      .attr("width", genWidth)
      .attr("height", 80);

    var margin = { right: 50, left: 50 };
    var width2 = +svg2.attr("width") - margin.left - margin.right;
    var height2 = +svg2.attr("height");

    var slider = svg2.append("g")
      .attr("class", "slider")
      .attr("transform", "translate(" + margin.left + "," + height2 / 2 + ")");

    slider.append("line")
      .attr("class", "track")
      .attr("x1", x.range()[0])
      .attr("x2", x.range()[1])
      .select(function () { return this.parentNode.appendChild(this.cloneNode(true)); })
      .attr("class", "track-inset")
      .select(function () { return this.parentNode.appendChild(this.cloneNode(true)); })
      .attr("class", "track-overlay")
      .call(d3.drag()
        .on("start.interrupt", function () { slider.interrupt(); })
        .on("start drag", function () {
          move(x.invert(d3.event.x));
        })
      );

    // label ticks
    slider.insert("g", ".track-overlay")
      .attr("class", "ticks")
      .attr("transform", "translate(0," + 18 + ")")
      .selectAll("text")
      .data(x.ticks(10))
      .enter()
      .append("text")
      .attr("x", x)
      .attr("y", 10)
      .attr("text-anchor", "middle")
      .text(function (d) { return d; });

    var handle = slider.insert("circle", ".track-overlay")
      .attr("class", "handle")
      .attr("r", 11);

    var cutoff = 0.0;

    function move(h) {
      cutoff = h;
      handle.attr("cx", x(cutoff));
      draw({ ...data, links: data.links.filter((link) => link.count > h) });
    }
  }

  return {
    remove,
    draw,
    createSlider,
  };
}