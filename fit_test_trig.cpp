#include "Eigen/Sparse"
#include "mesh_op.hpp"

int main()
{
    Mesh m1;

    if (m1.readmesh("magnet.msh"))
	    {
		    std::cout << "Mesh read: OK!" << std::endl;
	}
	else
    {
		std::cout << "Cannot open mesh file!" << std::endl;
        return 1;
    }
    std::vector<int> combo_list;
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
    for(int i = 0; i <m1.fg.size(); i++)// calculate face norm(orientation) & detect boundary type
    {
        m1.calc_face_norm(i);
        if(m1.fg[i].cells.size() == 2)
            if(m1.cg[m1.fg[i].cells[0]].tag == m1.cg[m1.fg[i].cells[1]].tag)
                m1.fg[i].boundary_type = 0;
            else
                m1.fg[i].boundary_type = 2;
        else
            m1.fg[i].boundary_type = 1;
    }    
    // apply loads to edges 
    for(int i = 0; i<m1.cg.size(); i++)
    {
        if(m1.cg[i].tag == 2)
        {
            for(int j = 0; j<m1.cg[i].faces.size(); j++)
                for(int k = 0; k<m1.fg[m1.cg[i].faces[j]].edges.size(); k++)
                {
                    int node1 = m1.eg[m1.fg[m1.cg[i].faces[j]].edges[k]].bnode;
                    int node2 = m1.eg[m1.fg[m1.cg[i].faces[j]].edges[k]].enode;
                    std::vector<double> edge_dir;
                    std::vector<double> load_dir={0.0,0.0,1.0}; //this is the value and direction of current density
                    for(int q = 0; q<3; q++)
                        edge_dir.push_back(m1.g[node2].glob_loc[q]-m1.g[node1].glob_loc[q]);
                    vec_norm(edge_dir);
                    sc_prod(edge_dir, load_dir, m1.eg[m1.fg[m1.cg[i].faces[j]].edges[k]].Load);
                }
        }
    }
// create a range of non-boundary edges
    std::vector<int> non_boundary_edges;
    for(int i = 0; i<m1.eg.size();i++)
    {
        for(int j=0; j<m1.eg[i].faces.size(); j++)
        {
            if(m1.fg[m1.eg[i].faces[j]].boundary_type != 1)
            {
                non_boundary_edges.push_back(i);
                m1.eg[i].range_position = non_boundary_edges.size()-1;
                break;
            }
        }
    }

    // just for test  will be removed at the end
    for(int i = 0; i <m1.fg.size(); i++)
    {
        for(int j=0; j<m1.fg[i].edges.size(); j++)
        {
            int node_num1 = m1.eg[m1.fg[i].edges[j]].bnode;
            int node_num2 = m1.eg[m1.fg[i].edges[j]].enode;
            m1.g[node_num1].phi += m1.fg[i].boundary_type;
            m1.g[node_num2].phi += m1.fg[i].boundary_type;
        }
    }
    for (int i = 0; i < m1.g.size(); i++)
    {
        for (int j = 0; j < m1.g[i].edges.size(); j++)
        {
            std::vector<double> edir;
            double value = m1.eg[m1.g[i].edges[j]].Load;
            m1.global_edge_dir(m1.g[i].edges[j], edir);
            m1.g[i].A[0] += edir[0]*value;
            m1.g[i].A[1] += edir[1]*value; 
            m1.g[i].A[2] += edir[2]*value; 
        }
        
    }

   // assemble matrix
    std::cout<<"start matrix assembly..."<<std::endl;

    typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
    typedef Eigen::Triplet<double> T;
    std::vector<T> coefficients;

    std::vector<int> tmp_edge_array;
// run over non_boundary edges only
    for(int i_nd=0; i_nd < non_boundary_edges.size();i_nd++)
    {
        int i = non_boundary_edges[i_nd]; //global index of edge
        double S = 0;
        std::vector<double> face_contrib;
        for(int j = 0; j < m1.eg[i].faces.size(); j++) // all faces of the edge
        {
            int face_number = m1.eg[i].faces[j];
            double dS; // part of the sectorial area
            double dl; // contribution to multiply all edges of the face
            face_contribution(i, face_number, dl, dS);
            S += dS;
            face_contrib.push(dl);
        }
        for(int j = 0; j < m1.eg[i].faces.size(); j++) // all faces of the edge
        {
            int face_number = m1.eg[i].faces[j];
            for(int k = 0; k < m.fg[face_number].edges.size(); k++)
            {
                int edge_number = m.fg[face_number].edges[k];

            }
            // determine edge orientation in the current face
            while (tmp_edge_array.size()>0)
            {
                for(int k = 0; k < tmp_edge_array.size(); k++)
                {
                    if(m.eg[tmp_edge_array[k]].bnode == next_node || m.eg[tmp_edge_array[k]].enode == next_node)
                    {
                        if(m.eg[tmp_edge_array[k]].bnode == next_node)
                        {
                            m.eg[tmp_edge_array[k]].sign = 1;
                            next_node = m.eg[tmp_edge_array[k]].enode;
                        }
                        else
                        {
                            m.eg[tmp_edge_array[k]].sign = -1;
                            next_node = m.eg[tmp_edge_array[k]].bnode;
                        }
                        tmp_edge_array.erase(tmp_edge_array.begin() + k); 
                        break;
                    }
                }
            }

            for(int k = 0; k < m.fg[m.eg[i].faces[j]].edges.size(); k++)
                {
                    int edge_index = m.fg[m.eg[i].faces[j]].edges[k];
                    if(edge_index != i)
                    {     
                        double coeff_value = m.eg[edge_index].sign/m.fg[m.eg[i].faces[j]].mu;
                        coefficients.push_back(T(i_nd,m.eg[edge_index].set_id, coeff_value));
                    }
                }
        }
    }
    SpMat Anb(non_boundary_edges.size(),non_boundary_edges.size());
    Anb.setFromTriplets(coefficients.begin(), coefficients.end());



    //m1.print_faces();
   // m1.writemesh("magnet_out.vtk");


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