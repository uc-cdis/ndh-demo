VERBOSITY_LEVEL=1;

fprintf (stdout, "\nRunning Bayesian Graphical Model discrete graph functionality for virus virulence...\n");

ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"bayesgraph.ibf");
ExecuteAFile ("/home/ubuntu/hyphy/virulence/createGBM.ibf");

/* Read binary matrix: 1 indicates headers 
   row and sets names and num_nodes 
*/
edgesAVG = {};
cv = 10;
for (k = 1; k < cv+1; k=k+1){
   
   fprintf (stdout, "Cross-validation simulation: "+k+"\n");
   mat = import_data("cv"+k+".csv", 1);

   nodes = {};
   for (i = 0; i < num_nodes; i=i+1) {

       /* add_discrete_node:
          node_id, max_parents, sample_size, nlevels
       */
       nodes[Abs(nodes)] = add_discrete_node(names[i], 5, 0, 2);
   }

   BayesianGraphicalModel bgm = (nodes);
   attach_data("bgm", mat, 0, 0, 0);
   res = graph_MCMC ("bgm", 1000000, 100000, 1000, 1);

   /* Write individual edges and calculate average */
   for (row = 0; row < num_nodes-1; row = row+1)
   {
	      for (col = row+1; col < num_nodes; col = col+1)
	      {
            edge = names[row]+","+names[col];
            edgesAVG[edge] = edgesAVG[edge] + res[row*num_nodes+col][1]/cv + res[col*num_nodes+row][1]/cv;
	      }
   }

}

create_edgelist("edgelist.out", edgesAVG, num_nodes, nodes);
create_mcmc_graph("bgm.dot", 0.5, edgesAVG, nodes);
