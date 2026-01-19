    fvScalarMatrix ZEqn
    (
        (
            fvm::ddt(rho, Z)
          + fvm::div(phi, Z)
          - fvm::laplacian(turbulence->DZEff(), Z)
        )
    );

    ZEqn.relax();
    //ZEqn.solve(mesh.solver("Z"));
    ZEqn.solve("Z");

    turbulence->correctChi();
    turbulence->correctVarZ();

    combustion->correct();
