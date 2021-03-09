function genomeGraph() {

  const svg = d3.select("div#genome")
    .append("svg:svg");

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

    const genomeGroup = genomesEnter
      .append("g")
      .attr("class", "genome")
      .enter();
    
    genomeGroup
      .append("svg:text")
      .text(function (d) { return d.accno; });

    genomeGroup
      .append("g")
      .attr("class", "proteins")
    
    const proteinsGroup = genomeGroup.selectAll("g.proteins").selectAll("rect")
      .data(function(d, i) {
        return d.proteins[i];
      })
      .enter();

    proteinsGroup
      .append("rect")
    
  }

  return {
    remove,
    draw, 
  }

}