!> \file canopy_app.F90
!! \brief Main Canopy Application Program
!! \details This is the main application program that coordinates the canopy model
!! calculations by orchestrating input reading, model initialization, time stepping,
!! canopy calculations, and output writing.
!!
!! \author Patrick C. Campbell
!! \date June 2022
!!
!! \version
!! - Prototype: Patrick C. Campbell, 06/2022
!! - Revised: PCC (10/2022)
!! - 21 Aug 2023: Adding multiple timesteps (P.C. Campbell)
!! - Revised: 30 Nov 2023: Added supplementary canopy profile, file_canvars (P.C. Campbell)

!> \defgroup canopy_app Canopy Application Main Program
!! \brief Main program for the canopy model application
!! \{

!> \brief Main Canopy Application Program
!! \details This program coordinates the entire canopy model workflow including:
!! - Reading user options from namelist
!! - Allocating and initializing variables
!! - Setting up input/output data structures
!! - Time stepping through the simulation period
!! - Calling main canopy calculations
!! - Writing output files
!! - Cleanup and deallocation
program canopy_app
    use canopy_date_mod    ! main canopy date module
    use canopy_files_mod   ! main canopy input files
    use canopy_coord_mod   ! main canopy coordinates

#ifdef NETCDF
    use canopy_ncf_io_mod
#endif

    implicit none

    !> \brief Next time step in YYYY-MM-DD-HH:MM:SS.SSSS format
    CHARACTER(LEN=24)                 :: time_next
    !> \brief Current time step in YYYY-MM-DD-HH:MM:SS.SSSS format
    CHARACTER(LEN=24)                 :: time_now
    !> \brief Current time formatted for file naming
    CHARACTER(LEN=24)                 :: time_now_file

    !> \brief Loop counters for time and general iteration
    integer   :: n,nn
!    CHARACTER(LEN=24)                 :: nn_string
!> \defgroup error_messages Error and Status Messages
!! \brief Formatted message strings for program output
!! \{

    !> \brief Format string for time processing header message
    CHARACTER(LEN=256), PARAMETER :: f100 = "(//, 1x, 78('~'), &
    & /,  1x, '~~~ Processing canopy-app for time = ', a,  &
    & /,  1x, 78('~'), /)"

!> \}

!> \defgroup main_workflow Main Application Workflow
!! \brief Primary computational workflow of the canopy application
!! \{


!> \brief Read user configuration options from namelist file
!! \details Reads the namelist.canopy file to get user-specified options
!! for the canopy model run including file paths, time settings, and
!! model parameters.
    call canopy_readnml

!> \brief Allocate memory for canopy model variables
!! \details Allocates all necessary arrays and data structures based on
!! the grid dimensions and model configuration read from the namelist.
    call canopy_alloc

!> \brief Initialize canopy model variables and parameters
!! \details Sets initial values for all canopy model variables, reads
!! biogenic emission parameters, and prepares the model state for calculations.
    call canopy_init

#ifdef NETCDF
    !> \brief Allocate NetCDF output data structures (NetCDF builds only)
    !! \details Allocates 2D/3D arrays for NetCDF output when the model
    !! is compiled with NetCDF support.
    call canopy_outncf_alloc

    !> \brief Initialize NetCDF output variables (NetCDF builds only)
    !! \details Initializes NetCDF output file structure and metadata
    !! when the model is compiled with NetCDF support.
    call canopy_outncf_init
#endif

!> \brief Main time stepping loop
!! \details Loops over the specified time period, reading input data,
!! performing canopy calculations, and writing output at each time step.
    time_now = time_start
    if(ntime.le.0) ntime=999999999 ! assign a large number
    timeloop: DO nn=1,ntime

        WRITE (*,f100) time_now
#ifdef NETCDF
        !> \brief Check and read NetCDF input files (NetCDF builds)
        !! \details Validates and reads meteorological and surface data from
        !! NetCDF input files, including main variables and canopy variables.
        call canopy_check_input(file_vars(nn),file_canvars(nn))
#else
        !> \brief Read text input files (non-NetCDF builds)
        !! \details Reads meteorological and surface data from text format
        !! input files when NetCDF support is not available.
        call canopy_read_txt(file_vars(nn),file_canvars(nn))
#endif

        !> \brief Perform main canopy model calculations
        !! \details Executes the core canopy model algorithms including:
        !! - Radiation calculations
        !! - Photosynthesis and emission calculations
        !! - Canopy profile computations
        !! - Biogenic emission estimates
        !! \param nn Current time step index
        call canopy_calcs(nn)

        !> \brief Format current time for output file naming
        !! \details Converts time string format by replacing colons with dashes
        !! to create valid filenames for output files.
        time_now_file=time_now
        do n = 1, len(time_now_file)
            if (time_now_file(n:n) == ':') then
                time_now_file(n:n) = '-'
            end if
        end do

        !> \brief Write text format output file
        !! \details Writes canopy model results to text format output file
        !! with timestamped filename.
        call canopy_write_txt((trim(file_out(1)) // '_' // trim(time_now_file)),time_now)

!        if (nn.lt.10) then
!            call canopy_write_txt((trim(file_out(1)) // '_t00' // ADJUSTL(nn_string)), &
!                time_now)
!        else if (nn.ge.10.and.nn.lt.100) then
!            call canopy_write_txt((trim(file_out(1)) // '_t0' // ADJUSTL(nn_string)), &
!                time_now)
!        else
!            call canopy_write_txt((trim(file_out(1)) // '_t' // ADJUSTL(nn_string)), &
!                time_now)
!        end if

#ifdef NETCDF
        !> \brief Write NetCDF format output file (NetCDF builds only)
        !! \details Writes 2D canopy model results to NetCDF format output file
        !! when NetCDF input files are used and NetCDF support is available.
        call canopy_write_ncf(trim(file_out(1)) // '_' // ADJUSTL(time_now_file))
#endif

        !> \brief Update to next time step
        !! \details Calculates the next time step based on the current time
        !! and the specified time interval.
        CALL geth_newdate (time_next, time_now, time_intvl)
        time_now = time_next

    ENDDO timeloop

!> \brief Deallocate memory for canopy model variables
!! \details Frees all allocated arrays and data structures to clean up
!! memory usage before program termination.
    call canopy_dealloc

!> \}

!> \}

    WRITE (*,'(//, a)') 'Canopy-App Finished Normally'

end program canopy_app
