function Graph() {
    var vertices = []; //存储图中所有的顶点名字
    var adjList = new Dictionary();//用之前的一个字典来存储邻接表
    this.addVertex = function(v){ //添加顶点
        vertices.push(v);
        adjList.set(v, []); //顶点为键，字典值为空数组
    };
    this.addEdge = function(v, w){ //添加边
        adjList.get(v).push(w); //基于有向图
        adjList.get(w).push(v); //基于无向图
    };
    this.toString = function(){
        var s = '';
        for (var i=0; i<vertices.length; i++){
            s += vertices[i] + ' -> ';
            var neighbors = adjList.get(vertices[i]);
            for (var j=0; j<neighbors.length; j++){
                s += neighbors[j] + ' ';
            }
            s += '\n';
        }
        return s;
    };
    var initializeColor = function(){
        var color = [];
        for (var i=0; i<vertices.length; i++){
            color[vertices[i]] = 'white';
        }
        return color;
    };
}  


// Show path function
// Only show vertex animation first
d3.select("#{{ figure_id }}_options .animate_button").on("click", function() {
    svg.selectAll(".metabolite")
    .transition()
    .style("fill", "#eb1818");
  
  
    });

//测试
var graph = new Graph();
var myVertices = ['A','B','C','D','E','F','G','H','I'];
for (var i=0; i<myVertices.length; i++){
    graph.addVertex(myVertices[i]);
}
graph.addEdge('A', 'B');
graph.addEdge('A', 'C');
graph.addEdge('A', 'D');
graph.addEdge('C', 'D');
graph.addEdge('C', 'G');
graph.addEdge('D', 'G');
graph.addEdge('D', 'H');
graph.addEdge('B', 'E');
graph.addEdge('B', 'F');
graph.addEdge('E', 'I');
console.log(graph.toString());