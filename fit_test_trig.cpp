#include "Eigen/Sparse"
#include "mesh_op.hpp"

int main()
{
    Mesh m1;
    std::vector<int> combo_list;

    if (m1.readmesh("test.msh"))
	    {
		    std::cout << "Mesh read: OK!" << std::endl;
	}
	else
		std::cout << "Cannot open mesh file!" << std::endl;


    make_combination_list(4,3,combo_list);

//try to make edges and faces 
    for (int i = 0; i < m1.cg.size(); i++)  
    {
        for (int j = 0; j<combo_list.size(); j++)
        {
            std::vector<int> bin_rep;
            std::vector<int> face_nodes;
            std::vector<int> face_edges;
            bool new_edge = false;
            int2bin(combo_list[j], bin_rep); // what nodes to take for the face as bin vector
            select_nodes(m1.cg[i].nodes, bin_rep, face_nodes); // the global numbers of these nodes

            int edge_number = m1.is_edge(face_nodes[0], face_nodes[1]); // test if thereis an edge
            if(edge_number == -1) // if there is no edge for 1-st pair of nodes
            {
                m1.eg.push_back(Edge(face_nodes[0],face_nodes[1])); //add edge to the list
                edge_number = m1.eg.size()-1; //update edge number
                m1.g[face_nodes[0]].edges.push_back(edge_number);//add edge to node
                m1.g[face_nodes[1]].edges.push_back(edge_number);
                new_edge=true;
            }
            face_edges.push_back(edge_number); // add edge number to list of edges of the current face
            edge_number = m1.is_edge(face_nodes[1], face_nodes[2]);
            if(edge_number == -1) // if there is no edge for 2-nd pair of nodes
            {
                m1.eg.push_back(Edge(face_nodes[1],face_nodes[2])); //add edge to the list
                edge_number = m1.eg.size()-1;
                m1.g[face_nodes[1]].edges.push_back(edge_number);
                m1.g[face_nodes[2]].edges.push_back(edge_number);
                new_edge=true;
            }
            face_edges.push_back(edge_number);
            edge_number = m1.is_edge(face_nodes[2], face_nodes[0]);
            if(edge_number == -1) // if there is no edge for 3-d pair of nodes
            {
                m1.eg.push_back(Edge(face_nodes[2],face_nodes[0])); //add edge to the list
                edge_number = m1.eg.size()-1;
                m1.g[face_nodes[2]].edges.push_back(edge_number);
                m1.g[face_nodes[0]].edges.push_back(edge_number);
                new_edge=true;
            }
            face_edges.push_back(edge_number);
            int face_number = m1.is_face(face_edges);
            if(face_number == -1) // if there is no face for this set of edges
            {
                m1.fg.push_back(Face(face_edges)); // add new face to mesh
                face_number = m1.fg.size()-1;
                for(int k = 0; k < face_edges.size(); k++)
                    m1.eg[face_edges[k]].faces.push_back(face_number);
            }
            // get face number to add to the cell
            //   std::cout<<"face edges "<<face_edges[0]<<" "<<face_edges[1]<<" "<<face_edges[2]<<std::endl;
            m1.cg[i].faces.push_back(face_number);
            m1.fg[face_number].cells.push_back(i);
        }   
    }
    for(int i = 0; i <m1.fg.size(); i++)
        m1.calc_face_norm(i);
    
    m1.print_faces();
    m1.writemesh("test_out.vtk");


 /* solving
    std::cout << "start BiCGSTAB ..." << std::endl;
 
    Eigen::BiCGSTAB<SpMat> solver;
    solver.compute(Anb);

    Eigen::VectorXd x = solver.solve(b);
    std::cout << "#iterations:     " << solver.iterations() << std::endl;
    std::cout << "estimated error: " << solver.error()      << std::endl;
    for (int i = 0; i < non_boundary_edges.size(); i++)
    {
        m.eg[non_boundary_edges[i]].Val = x(i);
    }
    
    //write vector potential in file
    std::ofstream MyFile("Afile.txt");

    for (int i = 0; i < m.g.size(); i++)
    {
        double Ax = 0.0;
        double Ay = 0.0;
        double Az = 0.0;

        for(int j = 0; j< m.g[i].edges.size();j++)
        {
            int edge_index = m.g[i].edges[j];
            double x0 = m.g[m.eg[edge_index].bnode].glob_loc[0];
            double y0 = m.g[m.eg[edge_index].bnode].glob_loc[1];
            double z0 = m.g[m.eg[edge_index].bnode].glob_loc[2];
            double x1 = m.g[m.eg[edge_index].enode].glob_loc[0];
            double y1 = m.g[m.eg[edge_index].enode].glob_loc[1];
            double z1 = m.g[m.eg[edge_index].enode].glob_loc[2];
            double L = std::sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+(z1-z0)*(z1-z0)); 
            Ax += std::abs(m.eg[edge_index].Val)*(x1-x0)/L;
            Ay += std::abs(m.eg[edge_index].Val)*(y1-y0)/L;
            Az += std::abs(m.eg[edge_index].Val)*(z1-z0)/L;
        }
        MyFile<<m.g[i].glob_loc[0]<<","<<m.g[i].glob_loc[1]<<","<<m.g[i].glob_loc[2]
        <<","<< Ax <<","<< Ay <<","<< Az <<std::endl; 
    }
    // Close the file
    MyFile.close();
    std::cout<<"Done!"<<std::endl;
*/
    return 0;
}