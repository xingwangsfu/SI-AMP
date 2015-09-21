# SI-AMP
SI-AMP algorithm

The code implements the proposed algorithm GENP-AMP in the paper “Approximate Message Passing-based Compressed Sensing Reconstruction with Generalized Elastic Net Prior” (http://arxiv.org/abs/1311.0576) and reproduces the experimental parts.

Requirements
1.	Mathworks MATLAB release 2009b or later
2.	The CVX software package, (available at http://cvxr.com/cvx/ )
3.	The GAMP software package, (available from from Sourceforge at http://sourceforge.net/projects/gampmatlab/files/), installed and included in MATLAB's path.

All related software packages are included in the document. If the readers want to get the latest version, please download from the webpages above.

Installation instructions:
1.	Install the CVX package from http://web.cvxr.com/cvx/doc/  and include the folders in MATLAB’s path.
2.	Install the GAMP software packages from http://sourceforge.net/projects/gampmatlab/files/  and include the folders in MATLAB’s path
3.	The multi-view dataset can be downloaded from http://www.fujii.nuee.nagoya-u.ac.jp/multiview-data/ .



How to use SI-OWLQN:
Save the data into .mtx format ,e.g., 
[ err1 ] = mmwrite( 'A_matrix.mtx',A); % %
[ err2 ] = mmwrite( 'y_matrix.mtx',Y); %
[ err3 ] = mmwrite( 'SI_matrix.mtx',x_SI); %
[ err4 ] = mmwrite('x_matrix.mtx',x);
The mmwrite and mmread programme have been included in sub-folder "solver".

cd to comparison algorithms\SI-OWLQN

The input format should be like this:
SI-OWLQN  x_matrix.mtx  A_matrix.mtx  y_matrix.mtx  SI_matrix.mtx  lambda  output.mtx  –l2weight  tau

Then mmread function is applied to read output.mtx into matlab and the MSE can be calculated.



