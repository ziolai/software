/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2017 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    laplacianFoam

Group
    grpBasicSolvers

Description
    Laplace equation solver for a scalar quantity.

    \heading Solver details
    The solver is applicable to, e.g. for thermal diffusion in a solid.  The
    equation is given by:

    \f[
        \ddt{T}  = \div \left( D_T \grad T \right)
    \f]

    Where:
    \vartable
        T     | Scalar field which is solved for, e.g. temperature
        D_T   | Diffusion coefficient
    \endvartable

    \heading Required fields
    \plaintable
        T     | Scalar field which is solved for, e.g. temperature
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
void ldu2csr(const Foam::lduMatrix & matrix, std::vector<int> & c_idx, std::vector<int> & r_idx,
	     std::vector<scalar> & vals); 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Laplace equation solver for a scalar quantity."
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating temperature distribution\n" << endl;

    fvScalarMatrix TEqn
    (
        - fvm::laplacian(DT, T)
           ==
          fvOptions(T)
    );

    fvOptions.constrain(TEqn);
    TEqn.solve();
    fvOptions.correct(T);

    const scalarField& diag = TEqn.diag();
    const scalarField& source = TEqn.source();
    const labelUList& neighbor = TEqn.lduAddr().lowerAddr();
    const labelUList& owner = TEqn.lduAddr().upperAddr();
    const scalarField& lower = TEqn.lower();
    const scalarField& upper = TEqn.upper();

    // if true, print LDU information on the matrix  
    if (1){
      Info<<" === Start LDU Information === "<< endl; 
      // print LDU diagonal matrix entries
      for(label i=0; i<diag.size(); i++)
      {
        Info<< "row index " << i << "\t" << " diag entry = "<< diag[i] << endl;
      }

      // print LDU off-diagonal entries
      for(label f=0; f<upper.size(); f++)
      {
        Info<< owner[f]    << "\t" << neighbor[f] << "\t" << lower[f] << endl;
        Info<< neighbor[f] << "\t" << owner[f]    << "\t" << upper[f] << endl;
      }
    
      // print right-hand side vector entries
      for(label i=0; i<diag.size(); i++)
      {
        Info<< "row index " << i << "\t" << " rhs entry = "<< source[i] << endl;
      }
      Info<<" === End LDU Information === "<< endl; 
    }
    
    // declare vectors to hold CSR matrix format
    std::vector< int >     rows;  // integer vector to hold row indices 
    std::vector< int >     cols;  // integer vector to hold column indices
    std::vector< scalar >  vals;  // real vector to hold matrix entries
                                  // scalar is here a type defined by Openfoam
    
    // allocate memory in vectors to hold CSR matrix format
    // note that the rows vector is allocated to problem size plus 2 as work memory 
    rows.reserve(diag.size()+2);
    cols.reserve(diag.size() + lower.size() + upper.size()); 
    vals.reserve(diag.size() + lower.size() + upper.size());
    // initialize rows vector to one: this is essential to proper working of ldu2csr() function 
    for(int i=0; i<diag.size()+2; i++) {
      rows[i] = 1;
    } ;

    // convert from LDU to CSR 
    ldu2csr(TEqn, cols, rows, vals); // mind the counter intuitive order of the arguments here  

    // if true, print LDU information on the matrix  
    if (1){
      Info<<" === Start CSR Information === "<< endl; 
    
      // print CSR cols indices
      for(int i=0; i<diag.size() + lower.size() + upper.size(); i++) {
        Info<< "  i = " << i  << "  cols[i] = " << cols[i] << endl;
      };
    
      // print CSR rows indices
      for(int i=0; i<=diag.size(); i++) {
        Info<< "  i = " << i  << "  rows[i] = " << rows[i] << endl;
      };
    
      // print CSR vals values 
      for(int i=0; i<diag.size() + lower.size() + upper.size(); i++) {
        Info<< "  i = " << i  << "  vals[i] = " << vals[i] << endl;
      };
       Info<<" === End CSR Information === "<< endl; 
    }
   

    #include "write.H"

    runTime.printExecutionTime(Info);
    

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
// Code provided  by BinWang0213
// at https://gist.github.com/BinWang0213/d5d4c9becf5d637cd4b8e25ac014cc48
#ifndef CSR_CONVERT_H
#define CSR_CONVERT_H

#include <vector>
#include "lduMatrix.H"

//
//	Conversion from OpenFOAM lduMatrix to classic Compressed Sparse Row matrix.
//
void ldu2csr( 
							const Foam::lduMatrix           & matrix, 
										 std::vector< int >       & c_idx, 
										 std::vector< int >       & r_idx,
										 std::vector< scalar > & vals
						) 
{
        int n_rows = matrix.diag().size();

	Info<<"  inside the ldu2csr function \n";
        Info<<"  ldu2csr():: diag.size()  = "<<  matrix.diag().size() << endl;
        Info<<"  ldu2csr():: lower.size() = "<<  matrix.lower().size() << endl;
        Info<<"  ldu2csr():: upper.size() = "<<  matrix.upper().size() << endl;

	//
	//	Calculate each row size. Sizes are shifted, because in the next part
	//	array r_idx is modified twice.
	//
	for (int i=0; i < matrix.upper().size() ; i++) {
 		int ri1 = matrix.lduAddr().lowerAddr()[i] ; 
		int ri2 = matrix.lduAddr().upperAddr()[i] ;
		r_idx[ri1+2] ++ ;
		r_idx[ri2+2] ++ ;
	} ;

	//
	//	Compute row offsets. Offsets are shifted by one positions, 
	//  because they are used as START positions while filling values.
	//
	r_idx[0] = 0 ;
	r_idx[1] = 0 ;
	for(int i=1; i<n_rows; i++) {
		r_idx[i+1] += r_idx[i] ;
	} ;
	
	//
	//	Fill in CSR matrix.
	//	Order below is important to keep column indices sorted.
	//

	// lower triangle
	for (int i=0; i < matrix.lower().size() ; i++) {
		int row    = matrix.lduAddr().upperAddr()[i] +1;
		int column = matrix.lduAddr().lowerAddr()[i] ;

		int idx = r_idx[row] ;
		vals[idx] = 0.; 
	        vals[idx] = matrix.lower()[i] ;
	        c_idx[idx] = column ;
		r_idx[row]++ ;
	} ;
	
	// diagonal
	for (int i=0; i<matrix.diag().size(); i++) {
		int idx = r_idx[i+1] ;
		vals[idx] = matrix.diag()[i] ;
		c_idx[idx] = i ; // i is row and column index
		r_idx[i+1]++ ;
	} ;

	// upper triangle
	for (int i=0; i < matrix.upper().size() ; i++) {
		int row    = matrix.lduAddr().lowerAddr()[i] +1;
		int column = matrix.lduAddr().upperAddr()[i] ;

		int idx = r_idx[row] ;
		vals[idx] = matrix.upper()[i] ;
		c_idx[idx] = column ;
		r_idx[row]++ ;
	} ;

} ;


#endif

// ************************************************************************* //


