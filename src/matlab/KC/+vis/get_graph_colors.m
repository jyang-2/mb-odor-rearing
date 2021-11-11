function node_coloring = get_graph_colors(FF2)
%GET_ Summary of this function goes here
%   Detailed explanation goes here

G = graph(FF2);
node_coloring = GraphColoringJohnson(G);

end

