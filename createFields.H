    int numPatches = mesh.boundaryMesh().size();
    wordList boundaryTypes(numPatches, "zeroGradient"); // default bc

    Info << "Create secondary (NO_READ) fields: " << endl;

    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "phi",
            dimArea*dimVelocity,
            scalar(0)
        )
    );

    volScalarField D
    (
        IOobject
        (
            "D",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "D", 
            dimensionSet(0,2,-1,0,0,0,0),
            scalar(1)
        ),
        boundaryTypes
    );

    volScalarField R
    (
        IOobject
        (
            "R",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "R", 
            dimensionSet(0,0,-1,0,0,0,0),
            scalar(1)
        ),
        boundaryTypes
    );

    volScalarField dR
    (
        IOobject
        (
            "dR",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "dR",
            dimensionSet(0,0,-1,0,0,0,0),
            scalar(1)
        ),
        boundaryTypes
    );

    List< autoPtr<volScalarField> > RList(nR);
    forAll(RList, n)
    {
        RList[n].reset
        (
            new volScalarField
            (
                IOobject
                (
                    "R" + name(n),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    "R" + name(n),
                    dimensionSet(0,0,-1,0,0,0,0),
                    scalar(1)
                ),
                boundaryTypes
            ) 
        );
    }

    List< autoPtr<volScalarField> > dRList(nR);
    forAll(dRList, n)
    {
        dRList[n].reset
        (
            new volScalarField
            (
                IOobject
                (
                    "dR" + name(n),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    "dR" + name(n),
                    dimensionSet(0,0,-1,0,0,0,0),
                    scalar(1)
                ),
                boundaryTypes
            ) 
        );
    }


    Info << "Reading primary fields: p, U, H, C" << endl;

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField H
    (
        IOobject
        (
            "H",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    List< autoPtr<volScalarField> > CList(nC);
    forAll(CList, n)
    {
        CList[n].reset
        (
            new volScalarField
            (
                IOobject
                (
                    "C" + name(n),
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            ) 
        );
    }

