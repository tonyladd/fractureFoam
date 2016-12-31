dimensionedScalar Uscale("Uscale", dimless, scalar(0));

Info << "Create velocity scale" << endl;
{
    volScalarField Htmp = H;            // Save current H

    if (runTime.value() != 0)
    {
        volScalarField H0               // Read H(t=0)
        (
            IOobject
            (
                "H",
                "0",
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        H = H0;
    }


//  Create Uscale 

    #include "solveU.H"


//  Restore current mineral fraction and velocity fields

    if (runTime.value() != 0)
    {
        H = Htmp;                   // Restore P
    }
}
