module canopy_models_mod

    implicit none

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE CANOPY_CANACC()

        !use canopy_modwr_mod, only: canacc_write
        !use canopy_modrd_mod, only: canacc_read

        integer :: run

        !write input 2D/3D canopy-app inputs to CANACC friendly input ./data files and general ./ctrl file containing file names
        !...canacc_write from module

        !run CANCACC inside canopy-app as separate submodule for each column
        call execute_command_line ('cd ./models/CANACC; ln -sf ./bin/canacc .; ./canacc', exitstat=run)
        print *, "Exit status of running the external model canacc was ", run

        !read each CANACC .out/ dat text files for each canopy column and variable
        !...canacc_read from module
        ! e.g., models/CANACC/out/'general_out'/canopy/tlsun.dat

        !store CANACC output variables into canopy-app arrays for each canopy column and variable, and pass out to main canopy-app
        !program...

    END SUBROUTINE CANOPY_CANACC

end module canopy_models_mod
