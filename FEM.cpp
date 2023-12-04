#include "inmost.h"
#include <stdio.h>


using namespace INMOST;
using namespace std;

/// Create three tags with different data types
/// and defined on different mesh elements
void main_create_tag(Mesh &m)
{
	Tag tagNodeCellNum;  // for each node stores number of surrounding cells
	Tag tagFaceArea;     // face areas
	Tag tagCellCentroid; // cell centroids

	// 1. Name: "NumberOfCells"
	// 2. Data type: integers
	// 3. Element type: mesh nodes
	// 4. Sparsity: not needed yet, set to NONE
	// 5. Tag size: 1
	tagNodeCellNum = m.CreateTag("NumberOfCells", DATA_INTEGER, NODE, NONE, 1);


	// 1. Name: "FaceArea"
	// 2. Data type: real numbers (double)
	// 3. Element type: mesh faces
	// 4. Sparsity: not needed yet, set to NONE
	// 5. Tag size: 1
	tagFaceArea = m.CreateTag("FaceArea", DATA_REAL, FACE, NONE, 1);

	// 1. Name: "CellCentroid"
	// 2. Data type: real numbers
	// 3. Element type: mesh cells
	// 4. Sparsity: not needed yet, set to NONE
	// 5. Tag size: 3
	tagCellCentroid = m.CreateTag("CellCentroid", DATA_REAL, CELL, NONE, 3);

	// Loops to fill tags

	// Node loop
	for(Mesh::iteratorNode inode = m.BeginNode(); inode != m.EndNode(); inode++){
		Node n = inode->getAsNode();
		int ncells = static_cast<int>(n.getCells().size());
		n.Integer(tagNodeCellNum) = ncells;
	}

	// Face loop
	for(Mesh::iteratorFace iface = m.BeginFace(); iface != m.EndFace(); iface++){
		Face f = iface->getAsFace();
		f.Real(tagFaceArea) = f.Area();
	}

	// Cell loop
	for(Mesh::iteratorCell icell = m.BeginCell(); icell != m.EndCell(); icell++){
		Cell c = icell->getAsCell();
		double x[3];
		c.Centroid(x);
		c.RealArray(tagCellCentroid)[0] = x[0];
		c.RealArray(tagCellCentroid)[1] = x[1];
		c.RealArray(tagCellCentroid)[2] = x[2];
	}
}
// Compute distance between two nodes
double nodeDist(const Node &n1, const Node &n2)
{
	double x1[3], x2[3];
	n1.Centroid(x1);
	n2.Centroid(x2);
	return sqrt(
				(  x1[0] - x2[0])*(x1[0] - x2[0])
				+ (x1[1] - x2[1])*(x1[1] - x2[1])
				+ (x1[2] - x2[2])*(x1[2] - x2[2])
				);
}

// Compute cell diameter
double cellDiam(const Cell &c)
{
	ElementArray<Node> nodes = c.getNodes();
	unsigned m = static_cast<unsigned>(nodes.size());
	double diam = 0.;
	for(unsigned i = 0; i < m; i++){
		for(unsigned j = 0; j < m; j++){
			diam = max(diam, nodeDist(nodes[i], nodes[j]));
		}
	}
	return diam;
}

// Compute and print mesh diameter
void main_mesh_diam(Mesh &m)
{
	double diam = 0.0;
	for(Mesh::iteratorCell icell = m.BeginCell(); icell != m.EndCell(); icell++){
		diam = max(diam, cellDiam(icell->self()));
	}
	cout << "Mesh diameter is " << diam << endl;
}
// Function to approximate
double g(double x, double y)
{
//return 0.0;
    double pi = 3.1415926535898;
    return sin(pi*x)*sin(pi*y);
}
double area(const Cell &c){
	return c.getFaces()[0].Area();
}
void center(const Cell & c, double * x){
	//c.Centroid(x);
    	ElementArray<Node> nodes = c.getNodes();
	x[0] = x[1] = x[2] = 0;
	for(unsigned i = 0; i < 3; i++){
		double coords[3];
		nodes[i].Centroid(coords);
		x[0] += coords[0] / 3;
		x[1] += coords[1] / 3;
	}

}

double basis_func(const Cell &c, const Node &n, double x_, double y_)
{
    ElementArray<Node> nodes = c.getNodes();
	unsigned n_ind = 0;
	double x[3];
	double y[3];
	for(unsigned i = 0; i < 3; i++){
		if(n == nodes[i])
			n_ind = i;
		double coords[3];
		nodes[i].Centroid(coords);
		x[i] = coords[0];
		y[i] = coords[1];
	}
	
	if(n_ind == 0){
		return ((x_   - x[2])*(y[1] - y[2]) - (x[1] - x[2])*(y_   - y[2])) /
			   ((x[0] - x[2])*(y[1] - y[2]) - (x[1] - x[2])*(y[0] - y[2]));
	}
	else if(n_ind == 1){
		return ((x_   - x[2])*(y[0] - y[2]) - (x[0] - x[2])*(y_   - y[2])) /
			   ((x[1] - x[2])*(y[0] - y[2]) - (x[0] - x[2])*(y[1] - y[2]));
	}
	else if(n_ind == 2){
		return ((x_   - x[0])*(y[1] - y[0]) - (x[1] - x[0])*(y_   - y[0])) /
			   ((x[2] - x[0])*(y[1] - y[0]) - (x[1] - x[0])*(y[2] - y[0]));
	}
	else{
		printf("Unexpected n_ind = %d\n", n_ind);
		exit(1);
	}
}
double dx_basis_func(const Cell &c, const Node &n, double x_, double y_)
{
    ElementArray<Node> nodes = c.getNodes();
	unsigned n_ind = 0;
	double x[3];
	double y[3];
	for(unsigned i = 0; i < 3; i++){
		if(n == nodes[i])
			n_ind = i;
		double coords[3];
		nodes[i].Centroid(coords);
		x[i] = coords[0];
		y[i] = coords[1];
	}
	cout << "dx_shit" << endl;
	for (int i = 0; i < 3; i++){
		cout << "(" << x[i] << ", " << y[i] << ")   ";
	}
	cout << endl;

	
	if(n_ind == 0){
		cout << (y[1] - y[2]) / ((x[0] - x[2])*(y[1] - y[2]) - (x[1] - x[2])*(y[0] - y[2])) << endl;
		return (y[1] - y[2]) / ((x[0] - x[2])*(y[1] - y[2]) - (x[1] - x[2])*(y[0] - y[2]));
	}
	else if(n_ind == 1){
		cout << (y[0] - y[2]) / ((x[1] - x[2])*(y[0] - y[2]) - (x[0] - x[2])*(y[1] - y[2])) << endl;
		return (y[0] - y[2]) / ((x[1] - x[2])*(y[0] - y[2]) - (x[0] - x[2])*(y[1] - y[2]));
	}
	else if(n_ind == 2){
		cout << (y[1] - y[0]) / ((x[2] - x[0])*(y[1] - y[0]) - (x[1] - x[0])*(y[2] - y[0])) << endl;
		return (y[1] - y[0]) / ((x[2] - x[0])*(y[1] - y[0]) - (x[1] - x[0])*(y[2] - y[0]));
	}
	else{
		printf("Unexpected n_ind = %d\n", n_ind);
		exit(1);
	}
}
double dy_basis_func(const Cell &c, const Node &n, double x_, double y_)
{
    ElementArray<Node> nodes = c.getNodes();
	unsigned n_ind = 0;
	double x[3];
	double y[3];
	for(unsigned i = 0; i < 3; i++){
		if(n == nodes[i])
			n_ind = i;
		double coords[3];
		nodes[i].Centroid(coords);
		x[i] = coords[0];
		y[i] = coords[1];
	}
	
	if(n_ind == 0){
		return (-(x[1] - x[2])) / ((x[0] - x[2])*(y[1] - y[2]) - (x[1] - x[2])*(y[0] - y[2]));
	}
	else if(n_ind == 1){
		return (-(x[0] - x[2])) / ((x[1] - x[2])*(y[0] - y[2]) - (x[0] - x[2])*(y[1] - y[2]));
	}
	else if(n_ind == 2){
		return (-(x[1] - x[0])) / ((x[2] - x[0])*(y[1] - y[0]) - (x[1] - x[0])*(y[2] - y[0]));
	}
	else{
		printf("Unexpected n_ind = %d\n", n_ind);
		exit(1);
	}
}
void main_linear_solver(Mesh &m)
{
	// Get number of nodes
	unsigned N = static_cast<unsigned>(m.NumberOfNodes());

	// Create sparse matrix, RHS vector and solution vector
	Sparse::Matrix A;
	Sparse::Vector b;
	Sparse::Vector sol;
	// Set their size - number of nodes
	A.SetInterval(0, N);
	b.SetInterval(0, N);
	sol.SetInterval(0, N);
	// Make A identity matrix
	// Make b: b_i = i
	for(Mesh::iteratorNode inode = m.BeginNode(); inode != m.EndNode(); inode++){
		Node n = inode->getAsNode();
		unsigned i = static_cast<unsigned>(n.LocalID());
		for (int j = 0; j < N; j++)
			A[i][j] = 0.0;
		b[i] = 0;
		ElementArray<Cell> arr = n.getCells();
		double areaSum = 0;
		unsigned m = static_cast<unsigned>(arr.size());
		for (int k = 0; k < m; k++){
			const Cell & c = arr[k];
			double cen[3];
			center(c, cen);
			ElementArray<Node> nodes = c.getNodes();
			unsigned m1 = static_cast<unsigned>(nodes.size());
			for (int k1 = 0; k1 < m1; k1++){
				const Node & nj = nodes[k1];
				unsigned j = static_cast<unsigned>(nj.LocalID());
				A[i][j] += area(c) * dx_basis_func(c, n, cen[0], cen[1]) * dx_basis_func(c, nj, cen[0], cen[1]);
				A[i][j] += area(c) * dy_basis_func(c, n, cen[0], cen[1]) * dy_basis_func(c, nj, cen[0], cen[1]);
			}
			areaSum += area(c);
			b[i] += area(c) * basis_func(c, n, cen[0], cen[1]) * g(cen[0], cen[1]);
		}
		for (int j = 0; j < N; j++){
			if (A[i][j] != 0)
				cout << "A[" << i << ", " << j << "] = " << A[i][j] << endl;
		}
		cout << "b[" << i << "] = " << b[i] << endl;
		cout << "area = " << areaSum << endl;
	}

	// Get solver
	// All inner INMOST solvers are BiCGStab
	// with different preconditioners, let's use ILU2
	Solver S(Solver::INNER_ILU2);
	S.SetParameter("absolute_tolerance", "1e-10");
	S.SetParameter("relative_tolerance", "1e-6");

	// Set matrix in the solver;
	// this also computes preconditioner
	S.SetMatrix(A);

	// Solve
	bool solved = S.Solve(b, sol);
	cout << "lin.it.: " << S.Iterations() << endl;
	if(!solved){
		cout << "Linear solver failure!" << endl;
		cout << "Reason: " << S.ReturnReason() << endl;
	}
	// Now we have solution in 'sol'
	// Create tag and save it here
	Tag tagSol = m.CreateTag("Sol", DATA_REAL, NODE, NONE, 1);
	for(Mesh::iteratorNode inode = m.BeginNode(); inode != m.EndNode(); inode++){
		Node n = inode->getAsNode();
		n.Real(tagSol) = sol[static_cast<unsigned>(n.LocalID())];
	}
}

int main(int argc, char ** argv)
{
	if( argc < 2 )
	{
		printf("Usage: %s mesh_file\n",argv[0]);
		return -1;
	}

	Mesh m;
	m.Load(argv[1]);
	main_create_tag(m);
	main_linear_solver(m);
	m.Save("res.vtk");
	cout << "end" << endl;
	return 0;
}
