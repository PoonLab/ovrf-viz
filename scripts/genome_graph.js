function genomeGraph() {
  const width = 800;
  const height = 800;

  const svg = d3.select("div#genome")
    .append("svg:svg")
    
    .attr("viewBox", `0 0 ${width} ${height}`);

  function remove() {
    svg.selectAll("g").selectAll('*').remove();
  }

  function draw(data) {
    
    // For each genome 
    const genomes = svg.selectAll("g.genome")
      .data(data, function(d){
        return d.accno;
      })
    
      // Create a group with class genome
    const genomesEnter = genomes
      .enter();

    const spacing = 2;  
    const genomeHeight = (height / data.length) - spacing;
    const genomeGroup = genomesEnter
      .append("g")
      .attr("class", "genome")
      .attr("transform", function(_d, i) {
        return `translate(0 ${i*(genomeHeight + spacing)})`;
      });
    
    // genomeGroup
    //   .append("svg:text")
    //   .attr("class", "textClass")
    //   .text(function (d) { return d.accno; });

    const proteinsGroup = genomeGroup
      .append("g")
      .attr("class", "proteins")
      
    
    const proteinsEnter = proteinsGroup.selectAll("rect")
      .data(function (d) {
        const genomeLength = d.proteins.reduce((maxEnd, p) => {
          return p.end > maxEnd ? p.end : maxEnd;
        }, -Infinity);

        return d.proteins.map((p) => ({
          ...p,
          startP: (p.start / genomeLength) * 100,
          endP: (p.end / genomeLength) * 100,
        }));
      })
      .enter();

    proteinsEnter
      .append("rect")
      .attr("x", function(d) { return `${d.startP}%` })
      .attr("y", 1)
      .attr("width", function(d) { return `${d.endP - d.startP}%` })
      .attr("height", genomeHeight)
      .attr("fill", function(d) { return d.color });
  }

  return {
    remove,
    draw, 
  }

}