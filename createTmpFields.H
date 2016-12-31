    volScalarField HTmp
    (
        IOobject
        (
            "HTmp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        H
    );

