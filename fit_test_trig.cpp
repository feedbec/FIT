#include "Eigen/Sparse"
#include "mesh_op.hpp"

//#define DIRECT

int main()
{
    Mesh m1;

    if (m1.readmesh("wire.msh"))
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
        for(int j=0; j<m1.cg[i].nodes.size(); j++)
            m1.g[m1.cg[i].nodes[j]].cells.push_back(i); // add cell number to node
    }
    int tp0 = 0, tp1=0, tp2 = 0;
    for(int i = 0; i <m1.fg.size(); i++)// calculate face norm(orientation) & detect boundary type
    {
        m1.calc_face_norm(i);
        double fdS;
        sc_prod(m1.fg[i].dS, m1.fg[i].dS, fdS);
        if(fdS< 1e-20)
            std::cout<<"very small face!"<<std::endl;
        
        if(m1.fg[i].cells.size() == 2){
            if(m1.cg[m1.fg[i].cells[0]].tag == m1.cg[m1.fg[i].cells[1]].tag)
            {
                m1.fg[i].boundary_type = 0;
                tp0++;
            }   
            else
            {
                m1.fg[i].boundary_type = 2;
                tp2++;
            }
        }
        else
        {
            m1.fg[i].boundary_type = 1;
            tp1++;
        }   
    }
    std::cout<<"outer f "<<tp1<<", inner f "<<tp2<<", other f "<<tp0<<std::endl;

    // apply material properties to cells
    double mu0 = 12.57e-7;
    for(int i = 0; i<m1.cg.size(); i++)
    {
        if(m1.cg[i].tag == 2)
        {
            m1.cg[i].mu = 1000.0*mu0;
        }
        else{
            m1.cg[i].mu = 1.0*mu0;
        }
    }
    

// create a range of non-boundary edges
    std::vector<int> non_boundary_edges;
    std::vector<int> boundary_edges;
    for(int i = 0; i<m1.eg.size();i++)
    {
        for(int j=0; j<m1.eg[i].faces.size(); j++)
        {
            if(m1.fg[m1.eg[i].faces[j]].boundary_type == 1)
            {
                boundary_edges.push_back(i);
                m1.eg[i].range_position = -2;
                break;
            }
        }
    }
    for(int i = 0; i<m1.eg.size();i++)
    {
        if(m1.eg[i].range_position == -1)
        {
            non_boundary_edges.push_back(i);
            m1.eg[i].range_position = non_boundary_edges.size()-1;
        }
    }

    // apply loads to edges analytically
    std::cout<<"start setting current"<<std::endl;
    for (int i = 0; i < m1.cg.size(); i++)
    {
        if(m1.cg[i].tag == 3){
            for (int j = 0; j < m1.cg[i].faces.size(); j++)
            {
                int face_number = m1.cg[i].faces[j];
                for (int k = 0; k < m1.fg[face_number].edges.size(); k++)
                {
                    int edge_number = m1.fg[face_number].edges[k];
                    std::vector<double> edge_center_vec;
                    m1.edge_center(edge_number, edge_center_vec);
                    edge_center_vec[0] = edge_center_vec[0] - 0.5;
                    edge_center_vec[1] = edge_center_vec[1] - 0.5;
                    std::vector<double> edge_dir;
                    double LD;
                    m1.edge_direction(edge_number, edge_dir);
                    std::vector<double> load_dir={0.0, 0.0, 1.0}; //this is the value and direction of current density
                    std::vector<double> load_dir1; 
                    vec_prod(load_dir, edge_center_vec, load_dir1);
                    vec_norm(load_dir1);
                    sc_prod(edge_dir, load_dir1, LD);
                    m1.eg[edge_number].Load = LD;
                }   
            }
        }
    }

   std::cout<<"Total edges "<<m1.eg.size()<<std::endl;
   std::cout<<"Non boundary edges "<<non_boundary_edges.size()<<std::endl;
   std::cout<<"Boundary edges "<<boundary_edges.size()<<std::endl;

   // assemble matrix
    std::cout<<"start matrix assembly..."<<std::endl;

    typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
    typedef Eigen::Triplet<double> T;
    std::vector<T> coefficients;
//    std::vector<int> tmp_edge_array;

// run over non_boundary edges only
    for(int i_nd=0; i_nd < non_boundary_edges.size();i_nd++)
    {
        all_edges pattern_edge_list;
        pattern_edge_list.clear();
        int i = non_boundary_edges[i_nd]; //global index of edge
        m1.rotrot(i, pattern_edge_list);
        m1.graddiv(i, pattern_edge_list);
        for (int j = 0; j < pattern_edge_list.edge_number.size(); j++)
        {
            coefficients.push_back(T(i_nd, m1.eg[pattern_edge_list.edge_number[j]].range_position,
            pattern_edge_list.edge_value[j]));
        }
    }

    SpMat Anb(non_boundary_edges.size(),non_boundary_edges.size());
    Anb.setFromTriplets(coefficients.begin(), coefficients.end());
 /*   for (int k=0; k<Anb.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(Anb,k); it; ++it)
        {
            std::cout<<it.value()<<" "<<it.row()<<" "<<it.col()<<std::endl;
    }
*/

    //set load vector
    Eigen::VectorXd b(non_boundary_edges.size()); 
    std::cout << "start loads..." << std::endl;
    for (int i = 0; i < non_boundary_edges.size(); i++)
    {
        b(i) = m1.eg[non_boundary_edges[i]].Load;
        //std::cout<<b(i)<<" ";
    }
    std::cout<<std::endl;

 // solving


#ifdef DIRECT
    std::cout << "start SparseLU ..." << std::endl;
     Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<int> >   solver;
    solver.analyzePattern(Anb); 
    // Compute the numerical factorization 
    solver.factorize(Anb); 
    //Use the factors to solve the linear system 
    Eigen::VectorXd x = solver.solve(b);  
#else
    std::cout << "start BiCGSTAB ..." << std::endl;
    Eigen::BiCGSTAB<SpMat> solver;
    solver.compute(Anb);
    solver.setMaxIterations(800);
    Eigen::VectorXd x = solver.solve(b);
    std::cout << "#iterations:     " << solver.iterations() << std::endl;
    std::cout << "estimated error: " << solver.error()      << std::endl;
#endif

    for (int i = 0; i < non_boundary_edges.size(); i++)
    {
        m1.eg[non_boundary_edges[i]].Val = x(i);
    }
    
    //interpolate vector potential to nodes
    for (int i = 0; i < m1.g.size(); i++)
        {
            double L = 0.0;
            for (int j = 0; j < m1.g[i].edges.size(); j++)
            {
                std::vector<double> edir;
                double dl;
                double value = m1.eg[m1.g[i].edges[j]].Val;
                m1.edge_direction(m1.g[i].edges[j], edir);
                sc_prod(edir,edir,dl);
                dl = 0.5*sqrt(dl);
                L += dl;
                vec_norm(edir);
                m1.g[i].A[0] += edir[0]*value*dl;
                m1.g[i].A[1] += edir[1]*value*dl; 
                m1.g[i].A[2] += edir[2]*value*dl; 
            }
            m1.g[i].A[0] /= L;
            m1.g[i].A[1] /= L; 
            m1.g[i].A[2] /= L; 
 
        }
    
    //interpolate current to nodes
    for (int i = 0; i < m1.g.size(); i++)
        {
            for (int j = 0; j < m1.g[i].edges.size(); j++)
            {
                std::vector<double> edir;
                double value = m1.eg[m1.g[i].edges[j]].Load;
                m1.edge_direction(m1.g[i].edges[j], edir);
                vec_norm(edir);
                m1.g[i].current[0] += edir[0]*value;
                m1.g[i].current[1] += edir[1]*value; 
                m1.g[i].current[2] += edir[2]*value; 
         } 
        }

    //interpolate B-field to cells
    for (int i = 0; i < m1.cg.size(); i++)
        {

            for (int j = 0; j < m1.cg[i].faces.size(); j++)
            {
                int face_number = m1.cg[i].faces[j];
                std::vector<double> rotA;
                std::vector<double> face_center_vec;
                m1.face_center(face_number, face_center_vec);
                double B = 0.0;
                for(int k=0;k<m1.fg[face_number].edges.size();k++)
                {
                    int edge_number = m1.fg[face_number].edges[k];
                    std::vector<double> edge_center_vec;
                    std::vector<double> face_edge;
                    std::vector<double> edir;
                    std::vector<double> orient;
                    m1.edge_center(edge_number, edge_center_vec);
                    for(int q=0;q<3;q++)
                        face_edge.push_back(edge_center_vec[q]-face_center_vec[q]);
                    m1.edge_direction(edge_number,edir);
                    vec_prod(face_edge,edir,orient);
                    vec_norm(orient);
                    double dS;
                    sc_prod(orient,m1.fg[face_number].dS,dS);
                    double dl;
                    sc_prod(edir,edir,dl);
                    dl = sqrt(dl);
                    B += dl*m1.eg[edge_number].Val/dS;
                }
                double dS;
                sc_prod(m1.fg[face_number].dS,m1.fg[face_number].dS,dS);
                dS = sqrt(dS);
                m1.cg[i].B[0] += B*m1.fg[face_number].dS[0]/dS;
                m1.cg[i].B[1] += B*m1.fg[face_number].dS[1]/dS;
                m1.cg[i].B[2] += B*m1.fg[face_number].dS[2]/dS;
            } 
        }


    //interpolate material properties to nodes
    for (int i = 0; i < m1.g.size(); i++)
        {
            m1.g[i].phi = (double) m1.cg[m1.g[i].cells[0]].tag;
        }

    std::cout << "Write results ..." << std::endl;

    m1.writemesh("wire_out_iter.vtk");
    std::cout<<"Done!"<<std::endl;
    return 0;
}