#include<vector>
#include<iostream>
#include<cmath>
#include <fstream>
#include <sstream>
#include <limits>

struct all_edges
{
    std::vector<int> edge_number;
    std::vector<double> edge_value;
    void add(int, double);
    void clear(){edge_number.clear(); edge_value.clear();}
};


struct Node
{
    double glob_loc[3];
    std::vector<int> edges;
    double phi;
    double A[3];
    double current[3];
    Node(double x, double y, double z)
    {glob_loc[0]=x;glob_loc[1]=y;glob_loc[2]=z; phi = 0.0; A[0]=A[1]=A[2]=0.0;
    current[0] = current[1] = current[2] = 0.0;};
};
typedef std::vector<Node> grid;

struct Edge
{
int bnode;
int enode;
double Val;
double Load;
std::vector<int> faces;
int range_position;
Edge(int n1, int n2):bnode(n1),enode(n2)
{Val = 0.0; Load = 0.0; range_position = -1;};
Edge(){Val = 0.0; Load = 0.0; range_position = -1;};
};
typedef std::vector<Edge> edge_grid;

struct Face
{
    std::vector<int> edges;
    std::vector<int> cells;
    std::vector<double> dS;
    int boundary_type; // 1-outer boundary, 0-inside block, 2-outer boundary
    Face(){};
    Face(const std::vector<int>&a):edges(a){};
};
typedef std::vector<Face> face_grid;

struct Cell
{
    double mu;
    int tag;
    std::vector<int> nodes;
    std::vector<int> faces;
    Cell(){};
    Cell(std::vector<int> & nodes0, int tag0):nodes(nodes0), tag(tag0){mu = 1.0;};
};
typedef std::vector<Cell> cell_grid;

struct Mesh
{
    grid g;
    edge_grid eg;
    face_grid fg;
    cell_grid cg;
    Mesh(){};
    void print_nodes();
    void print_edges();
    void print_faces();
    bool readmesh(std::string filename) ;
    bool writemesh(std::string filename);
    void common_edges(const std::vector<int>& node_list, std::vector<int>& edge_list);
    int is_face(const std::vector<int>& edge_list);
    int is_edge(int N1, int N2); 
    void calc_face_norm(int i);
    void edge_direction(int edge_num, std::vector<double> & v);
    void face_contribution(int edge_number,
        int face_number, double & face_contrib, double & sec_S);
    void edge_contribution(int edge_num, int face_num, double & dl);
    void face_center(int face_number, std::vector<double> & face_center_vec);
    void edge_center(int edge_number, std::vector<double> & edge_center_vec) ;
    void cell_center(int cell_number, std::vector<double> & cell_center_vec) ;
    void rotrot(int edge_number, all_edges & elist);
    void graddiv(int edge_number, all_edges & elist);
    double dual_face_area(int edge_number);
    double dual_cell_volume(int node_number);
};

std::vector<std::string> split(const std::string &s, char delim);
void int2bin(int N, std::vector<int> & bin_representation);
void select_nodes(const std::vector<int> & basic_set, const std::vector<int> & to_take,
std::vector<int> & result );
void vec_prod05(const std::vector<double> & v1, const std::vector<double> & v2,
std::vector<double> & dS);
void vec_prod(const std::vector<double> & v1, const std::vector<double> & v2,
std::vector<double> & dS);
void sc_prod(const std::vector<double> & v1, const std::vector<double> & v2, double & result);
void vec_norm(std::vector<double> & v2);

void Mesh::print_faces()
{
    for (int i = 0; i < fg.size(); i++)
    {
        std::cout<< "Face "<< i <<":"<<std::endl;
        std::cout<< "Edges ";
        for (int j = 0; j < fg[i].edges.size(); j++)
        {
            std::cout<<fg[i].edges[j]<<" ";
        }
        std::cout<<std::endl;
        std::cout<<"Cells ";
        for (int j = 0; j < fg[i].cells.size(); j++)
        {
            std::cout<<fg[i].cells[j]<<" ";
        }
        std::cout<<std::endl;

    }
    
}
void Mesh::print_edges()
{
    for (int i = 0; i < eg.size(); i++)
    {
        std::cout<<i<<" nodes: "<<eg[i].bnode<<" "<<eg[i].enode<<" faces: ";
        for (int j = 0; j < eg[i].faces.size(); j++)
        {
            std::cout<<eg[i].faces[j]<<" ";
        }
        std::cout<<std::endl;
    }
    
}

void Mesh::print_nodes()
{
    for (int i = 0; i < g.size(); i++)
    {
        std::cout<<i<<" "<<g[i].glob_loc[0]<<" "
        <<g[i].glob_loc[1]<<" "<<g[i].glob_loc[2]<<" ";
        for(int j=0; j<g[i].edges.size(); j++)
            std::cout<<g[i].edges[j]<<" ";
        std::cout<<std::endl;
    }
    
}

void Mesh::common_edges(const std::vector<int>& node_list, std::vector<int>& edge_list)
{   
    edge_list.clear();
    std::vector<int> pool;
    for(int i = 0; i< g[node_list[0]].edges.size(); i++)
        pool.push_back(g[node_list[0]].edges[i]);
    
    for(int i=1; i<node_list.size();i++)
    {
        int pool_size = pool.size();
        for(int j = 0; j <g[node_list[i]].edges.size(); j++)
        {
            bool in_pool = false;
            for(int k = 0; k< pool_size; k++)
            {
                if (g[node_list[i]].edges[j] == pool[k])
                {
                    edge_list.push_back(pool[k]);
                    in_pool = true;
                    break;
                }
            }
            if(!in_pool)
            {
                pool.push_back(g[node_list[i]].edges[j]);
            }

        }
    }
}

int Mesh::is_face(const std::vector<int>& edge_list)
{   
    for(int i = 0; i< eg[edge_list[0]].faces.size(); i++)
        for(int j = 0; j< eg[edge_list[1]].faces.size(); j++)
            if(eg[edge_list[0]].faces[i] == eg[edge_list[1]].faces[j])
            {
                return eg[edge_list[1]].faces[j];
            }
    return -1;
}

//generate list of subsets of K elements from total N
//combination list contains decimals that after conversion to bin 
//indicates which element to take
void make_combination_list(int N, int K, std::vector<int> & combination_list)
{
    combination_list.clear();
    std::vector<int> comb;
    int S = 0;
    for(int i = N-1; i >= 0; i--)
    {
        int S0 = 1;
        for(int j = 0; j<i; j++)
            S0 *=2;
        S += S0;
    }
    for(int i = 0; i<S;i++)
    {
        //convert to binary and calculate occupied positions
        int  nops = 0;
        int2bin(i, comb);
        for (int j = 0; j < comb.size(); j++)
        {
            nops += comb[j];
        }
        
        if(nops != K)
            continue;
        else
            combination_list.push_back(i);

    }
}

// decimal to binary converter bin_representation contains 0 and 1
void int2bin(int N, std::vector<int> & bin_representation)
{
    bin_representation.clear();
    int  div, val = N, rest;
    do{
        div = val / 2;
        rest = val % 2;
        bin_representation.push_back(rest);
        val = div; 
    }while (div > 0);  
}

void select_nodes(const std::vector<int> & basic_set, const std::vector<int> & to_take,
std::vector<int> & result )
{
    result.clear();
    for (int i = 0; i < to_take.size(); i++)
    {
        if(to_take[i])
            result.push_back(basic_set[i]);
    }
}

int Mesh::is_edge(int N1, int N2)
{
    for(int i = 0; i<g[N1].edges.size(); i++)
        for(int j = 0; j<g[N2].edges.size(); j++)
            if(g[N1].edges[i] == g[N2].edges[j])
                return g[N1].edges[i];
    return -1;
}

bool Mesh::writemesh(std::string filename)
{
	std::fstream newfile;

	newfile.open(filename, std::ios::out); //open a file to perform read operation using file object
	if (newfile.is_open())
	{
		newfile << "# vtk DataFile Version 1.0" << std::endl;
		newfile << "Unstructured Grid Example" << std::endl;
		newfile << "ASCII" << std::endl<< std::endl;
		newfile << "DATASET UNSTRUCTURED_GRID" << std::endl;
		newfile << "POINTS "<< g.size() <<" float" << std::endl;

		for (size_t i = 0; i < g.size(); i++)
			newfile << g[i].glob_loc[0] << " " << g[i].glob_loc[1] << " " << g[i].glob_loc[2] << std::endl;

		newfile << "CELLS " << cg.size() << " " << 5 * cg.size() << std::endl;
		for (size_t i = 0; i < cg.size(); i++)
		{
			newfile << 4 << " ";
			for (size_t j = 0; j < 4; j++)
				newfile << cg[i].nodes[j] <<" ";
			newfile << std::endl;
		}
		newfile << std::endl;
		newfile << "CELL_TYPES " << cg.size() << std::endl;
		for (size_t i = 0; i < cg.size(); i++)
			newfile << 10 << std::endl;
		newfile << std::endl;
		newfile << "POINT_DATA " << g.size() << std::endl;
		newfile << "SCALARS scalar_potential float" << std::endl;
		newfile << "LOOKUP_TABLE default" << std::endl;
		for (size_t i = 0; i < g.size(); i++)
			newfile << g[i].phi << std::endl;
		// print vector data
        newfile << std::endl;
		newfile << "VECTORS current float" << std::endl;
		for (size_t i = 0; i < g.size(); i++)
			newfile << g[i].current[0]<<" "<< g[i].current[1]<<" " << g[i].current[2] << std::endl;	

        newfile << std::endl;
		newfile << "VECTORS Vector_potential float" << std::endl;
		for (size_t i = 0; i < g.size(); i++)
			newfile << g[i].A[0]<<" "<< g[i].A[1]<<" " << g[i].A[2] << std::endl;	
        newfile.close();
		return true;
	}
	else
		return false;
}

//void Mesh::edge2node(const std::vector<double> & Edge_data, std::vector<vector<double>> & Node_data)


bool Mesh::readmesh(std::string filename) 
{
	bool node_read = false;
	bool elem_read = false;
	int node_max;
	int elem_max;
	int tag = 0;
	std::fstream newfile;

	newfile.open(filename, std::ios::in); //open a file to perform read operation using file object
	if (newfile.is_open()) {   //checking whether the file is open
		std::string tp;
		while (getline(newfile, tp)) { //read data from file object and put it into string.
			if (tp == "$EndNodes") // start node reading
				node_read = false;
			if (node_read)
			{
				std::vector<std::string> v = split(tp, ' ');
				if ((int)v.size() == 3)
				{
					g.push_back(Node(stod(v[0]), stod(v[1]), stod(v[2])));
				}
			}
			if (tp == "$Nodes")
			{
				node_read = true;
				getline(newfile, tp);
				std::vector<std::string> v = split(tp, ' ');
				node_max = stoi(v[3]);
			}

			if (tp == "$EndElements") // start element reading
				elem_read = false;
			if (elem_read)
			{
				std::vector<std::string> v = split(tp, ' ');
				if ((int)v.size() == 5)
				{
					std::vector<int> temp_set = { stoi(v[1]) - 1, stoi(v[2]) - 1, stoi(v[3]) - 1, stoi(v[4]) - 1 };
					std::sort(temp_set.begin(), temp_set.end());
					cg.push_back(Cell(temp_set,tag));
				}
				
				if ((int)v.size() == 4)
				{
					tag = stoi(v[1]);
				}


			}
			if (tp == "$Elements")
			{
				elem_read = true;
				getline(newfile, tp);
				std::vector<std::string> v = split(tp, ' ');
				elem_max = stoi(v[3]);
			}

		}
		newfile.close(); //close the file object.
		std::cout << "Total number of nodes: " << g.size() << " Expected: " << node_max << std::endl;
		std::cout << "Total number of elements: " << cg.size() << " Expected: " << elem_max << std::endl;

		return true;
	}
	else
		return false;
}

std::vector<std::string> split(const std::string &s, char delim) {
	std::vector<std::string> result;
	std::stringstream ss(s);
	std::string item;

	while (getline(ss, item, delim)) {
		result.push_back(item);
	}

	return result;
}

void Mesh::calc_face_norm(int i)
{
    std::vector<double> v1, v2;
    edge_direction(fg[i].edges[0],v1);
    edge_direction(fg[i].edges[1],v2);
    vec_prod05(v1,v2,fg[i].dS);    
}

void vec_prod(const std::vector<double> & v1, const std::vector<double> & v2,
std::vector<double> & dS)
{
    dS.clear();
    dS.push_back(v1[1]*v2[2]-v1[2]*v2[1]);
    dS.push_back(-v1[0]*v2[2]+v1[2]*v2[0]);
    dS.push_back(v1[0]*v2[1]-v1[1]*v2[0]);
}

void vec_prod05(const std::vector<double> & v1, const std::vector<double> & v2,
std::vector<double> & dS)
{
    dS.clear();
    dS.push_back(0.5*(v1[1]*v2[2]-v1[2]*v2[1]));
    dS.push_back(0.5*(-v1[0]*v2[2]+v1[2]*v2[0]));
    dS.push_back(0.5*(v1[0]*v2[1]-v1[1]*v2[0]));
}

void sc_prod(const std::vector<double> & v1, const std::vector<double> & v2, double & result)
{
    result = 0;
    for(int i = 0; i< v1.size(); i++)
        result += v1[i]*v2[i];
}

void vec_norm(std::vector<double> & v2)
{
    double norm = 0.0;
    const double TOL = std::numeric_limits< double >::min();

    sc_prod(v2,v2,norm);
    if(norm < TOL)
    {
        std::cout<<"zero norm"<<std::endl;
        return;
    }

    for (int i = 0; i < v2.size(); i++)
    {
        v2[i] = v2[i]/sqrt(norm);
    }   
}

void Mesh::face_center(int face_number, std::vector<double> & face_center_vec)
{
    double q1, q2;
    face_center_vec.clear();
    for(int i = 0; i<3; i++)
        face_center_vec.push_back(0.0);
    int N = fg[face_number].edges.size(); //number of edges
    for(int i = 0; i< N; i++)
    {
        std::vector<double> edge_center_vec;
        edge_center(fg[face_number].edges[i],edge_center_vec);
        for(int j=0; j<3;j++)
        {
            face_center_vec[j] += edge_center_vec[j]/N;
        }
    }
}

void Mesh::edge_center(int edge_number, std::vector<double> & edge_center_vec)
{
    double q1,q2;
    edge_center_vec.clear();
    for(int i = 0; i<3; i++)
    {
        q1 = g[eg[edge_number].bnode].glob_loc[i];
        q2 = g[eg[edge_number].enode].glob_loc[i];
        edge_center_vec.push_back(0.5*(q2+q1));
    }
}

void Mesh::edge_direction(int edge_number, std::vector<double> & edge_dir_vec)
{
    double q1,q2;
    edge_dir_vec.clear();
    for(int i = 0; i<3; i++)
    {
        q1 = g[eg[edge_number].bnode].glob_loc[i];
        q2 = g[eg[edge_number].enode].glob_loc[i];
        edge_dir_vec.push_back((q2-q1));
    }
}

void Mesh::cell_center(int cell_number, std::vector<double> & cell_center_vec)
{
    cell_center_vec.clear();

    for(int i = 0; i<3; i++)
    {
        double q = 0.0;
        for(int j=0; j<cg[cell_number].nodes.size();j++)
        {
            q += g[cg[cell_number].nodes[j]].glob_loc[i];
        }
        cell_center_vec.push_back((q/cg[cell_number].nodes.size()));
    }
}

void Mesh::face_contribution(int edge_number,
    int face_number, double & face_contrib, double & sec_S)
{
    if(fg[face_number].cells.size() != 2)
    {
        std::cout<<"Not 2 cells"<<std::endl;
        return;
    }
    std::vector<double> edir_vec;
    std::vector<double> fdir_vec;
    std::vector<double> edge_center_vec;
    std::vector<double> face_center_vec;
    std::vector<double> cell_center_vec;
    std::vector<double> edge_face;
    std::vector<double> cell_face;
    std::vector<double> sec_vec;

    face_center(face_number, face_center_vec);
    edge_center(edge_number, edge_center_vec);
    edge_direction(edge_number, edir_vec);

    //from edge center to face center
    for (int i = 0; i < 3; i++)
    {
        edge_face.push_back(face_center_vec[i] - edge_center_vec[i]);
    }
    vec_prod(edir_vec, edge_face, fdir_vec);
    vec_norm(fdir_vec);
    double dS; //sign and area of the face corresponding to the given edge rotation direction
    sc_prod(fdir_vec, fg[face_number].dS, dS);
    
    // calculate length of the circulation path in 1-st cell
    int cell_number = fg[face_number].cells[0];
    double dl;
    std::vector<double> cell_to_face;
    cell_center(cell_number, cell_center_vec); 
    // distance from cell center to face center
    for (int i = 0; i < 3; i++)
    {
        cell_face.push_back( face_center_vec[i]-cell_center_vec[i]);
    }
    sc_prod(cell_face, cell_face, dl);
    dl = sqrt(dl)/dS;
    //dl = dS/fabs(dS);
    double cont1 = dl/(cg[cell_number].mu);
    //calculate the sectorial area in 1-st cell
    vec_prod05(cell_face, edge_face, sec_vec);
    double sec_S1;
    sc_prod(sec_vec, sec_vec, sec_S1);
    sec_S1 = sqrt(sec_S1);

    //distance to 2-nd cell center
    cell_number = fg[face_number].cells[1];
    cell_center(cell_number, cell_center_vec); 
    // distance from cell center to face center
    cell_face.clear();
    for (int i = 0; i < 3; i++)
    {
        cell_face.push_back( face_center_vec[i]-cell_center_vec[i]);
    }
    sc_prod(cell_face,cell_face, dl);
    dl = sqrt(dl)/dS;
    //dl = dS/fabs(dS);
    double cont2 = dl/(cg[cell_number].mu);
    //calculate the sectorial area in 1-st cell
    vec_prod05(cell_face, edge_face, sec_vec);
    double sec_S2;
    sc_prod(sec_vec, sec_vec, sec_S2);
    sec_S2 = sqrt(sec_S2);
    face_contrib = cont1+cont2; 
    sec_S = sec_S1 + sec_S2;
}


void Mesh::edge_contribution(int edge_number, int face_number, double & dl)
{
    std::vector<double> edge_center_vec;
    std::vector<double> face_center_vec;
    std::vector<double> edir_vec;
    std::vector<double> face_edge;
    std::vector<double> orient;

    edge_center(edge_number, edge_center_vec);
    face_center(face_number, face_center_vec);
    edge_direction(edge_number, edir_vec);
    sc_prod(edir_vec, edir_vec, dl);
    dl = sqrt(dl);
    for (int i = 0; i < 3; i++)
    {
        face_edge.push_back(-edge_center_vec[i] + face_center_vec[i]);
    }
    vec_prod(edir_vec, face_edge, orient);
    vec_norm(orient);
    double signum;
    sc_prod(fg[face_number].dS, orient, signum);
    signum = signum/fabs(signum);
    dl = signum*dl; //
//    std::cout<<"EC done"<<std::endl;
}

void Mesh::rotrot(int edge_number, all_edges & elist)
{
    int i = edge_number;
    double S = 0.0; //total area around the edge
    std::vector<double> face_contrib;
    face_contrib.clear();
    for(int j = 0; j < eg[i].faces.size(); j++) // all faces of the edge
    {
        int face_number = eg[i].faces[j];
        double dS; // part of the sectorial area
        double dl; // contribution to multiply all edges of the face
        face_contribution(i, face_number, dl, dS);
        S += dS;
        face_contrib.push_back(dl);
    }
    for(int j = 0; j < face_contrib.size(); j++) // all faces of the edge
    {
        face_contrib[j] = face_contrib[j]/S;
    }

    for(int j = 0; j < eg[i].faces.size(); j++) // all faces of the edge
    {
        int face_number = eg[i].faces[j];
        for(int k = 0; k < fg[face_number].edges.size(); k++) // all edges of the face
        {
            int edge_number = fg[face_number].edges[k];
            double dl;
            edge_contribution(edge_number, face_number, dl);
            double coeff_value = dl * face_contrib[j];
            if(eg[edge_number].range_position >=0) 
                elist.add(edge_number, coeff_value);
        }
    }

}

double Mesh::dual_face_area(int edge_number)
{
    double S = 0.0; //total area around the edge
    face_contrib.clear();
    for(int j = 0; j < eg[edge_number].faces.size(); j++) // all faces of the edge
    {
        int face_number = eg[edge_number].faces[j];
        double dS; // part of the sectorial area
        double dl; // contribution to multiply all edges of the face
        face_contribution(edge_number, face_number, dl, dS);
        S += dS;
    }
    return S;
}

double Mesh::dual_cell_volume(int node_number)
{

}


void all_edges::add(int num, double val)
{
    for(int i = 0; i <edge_number.size(); i++)
    {
        if(num == edge_number[i])
        {
            edge_value[i] += val;
            return;
        }
    }
    edge_number.push_back(num);
    edge_value.push_back(val);
}

void Mesh::graddiv(int edge_number, all_edges & elist)
{
    int i = edge_number;
    double S = 0.0; //total area around the edge
    std::vector<double> face_contrib;
    face_contrib.clear();
    for(int j = 0; j < eg[i].faces.size(); j++) // all faces of the edge
    {
        int face_number = eg[i].faces[j];
        double dS; // part of the sectorial area
        double dl; // contribution to multiply all edges of the face
        face_contribution(i, face_number, dl, dS);
        S += dS;
        face_contrib.push_back(dl);
    }
    for(int j = 0; j < face_contrib.size(); j++) // all faces of the edge
    {
        face_contrib[j] = face_contrib[j]/S;
    }

    for(int j = 0; j < eg[i].faces.size(); j++) // all faces of the edge
    {
        int face_number = eg[i].faces[j];
        for(int k = 0; k < fg[face_number].edges.size(); k++) // all edges of the face
        {
            int edge_number = fg[face_number].edges[k];
            double dl;
            edge_contribution(edge_number, face_number, dl);
            double coeff_value = dl * face_contrib[j];
            if(eg[edge_number].range_position >=0) 
                elist.add(edge_number, coeff_value);
        }
    }

}

