

    Info << "***************************************************" << nl;
    Info << "Reading fractureFoam dictionary" << nl;
    IOdictionary fractureFoam
    (
        IOobject
        (
            "fractureFoamDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
                                                       
    int nC;
    if( !fractureFoam.readIfPresent<int>("nC", nC) ){
    SeriousErrorIn("main")
    << "There is no nC parameter in fractureFoam dictionary"
    << exit(FatalError);
    }
    int nR;
    if( !fractureFoam.readIfPresent<int>("nR", nR) ){
    SeriousErrorIn("main")
    << "There is no nR parameter in fractureFoam dictionary"
    << exit(FatalError);
    }

    bool QCON;
    if( !fractureFoam.readIfPresent<bool>("QCON", QCON) ){
    SeriousErrorIn("main")
    << "There is no QCON parameter in fractureFoam dictionary"
    << exit(FatalError);
    }
    bool DEBUG;
    if( !fractureFoam.readIfPresent<bool>("DEBUG", DEBUG) ){
    SeriousErrorIn("main")
    << "There is no DEBUG parameter in fractureFoam dictionary"
    << exit(FatalError);
    }

    double Qmax;
    if( !fractureFoam.readIfPresent<double>("Qmax", Qmax) ){
    SeriousErrorIn("main")
    << "There is no Qmax parameter in fractureFoam dictionary"
    << exit(FatalError);
    }

    Info << "FractureFoamDict:            "           << nl;
    Info << "# of concentrations :        " << nC     << nl;
    Info << "# of reactions :             " << nR     << nl;
    Info << "Constant Flow Rate:          " << QCON   << nl;
    Info << "Debugging:                   " << DEBUG  << nl;
    Info << "Max Flow Rate (QCON false):  " << Qmax   << nl;

    Info << "Reading fvSolution controls dictionary\n" << endl;
    const dictionary& controls = mesh.solutionDict().subDict("concentration");

    int nCycle;
    if( !controls.readIfPresent<int>("nCycle", nCycle) ){
    SeriousErrorIn("main")
    << "There is no nCycle parameter in fractureFoam dictionary"
    << exit(FatalError);
    }
    double tolerance;
    if( !controls.readIfPresent<double>("tolerance", tolerance) ){
    SeriousErrorIn("main")
    << "There is no tolerance parameter in fractureFoam dictionary"
    << exit(FatalError);
    }

    Info << "fvSolution.concentrationDict:"              << nl;
    Info << "# cycles:                    " << nCycle    << nl;
    Info << "tolerance:                   " << tolerance << nl << nl;

    Info << "Reading transportProperties\n" << nl;
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar D0 (transportProperties.lookup("D0"));
    dimensionedScalar k0 (transportProperties.lookup("k0"));

    Info << "Constant difffusion coefficient in fractures" << endl;
    Info << "D0 = U0*L/Pe: " << D0 << endl;
    Info << "k0 = Da*U0/L: " << k0 << endl;

    Info << "***************************************************" << nl
         << endl;

