//  Calculate R and dR for a species; include mineral conditions

{
    #include "reacModel.H"              // Get R & dR for all reactions

    dimensionedScalar Zero("", dimless/dimTime, scalar(0));
    R = Zero;
    dR = Zero;
    for (int r=0; r < nR; r++){
        R +=  RList[r]()*stoiTable[r][Case];
       dR += dRList[r]()*stoiTable[r][Case];
    }
     R.correctBoundaryConditions();
    dR.correctBoundaryConditions();
}
