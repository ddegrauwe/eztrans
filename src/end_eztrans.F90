subroutine cleanup(config)

use mpi

use config_mod, only : config_type

implicit none


! arguments
type(config_type), intent(inout) :: config

! deallocations
deallocate(config%nx_l)
deallocate(config%jxi_l)
deallocate(config%jxe_l)
deallocate(config%ny_l)
deallocate(config%jyi_l)
deallocate(config%jye_l)
deallocate(config%nfld_l)
deallocate(config%jfldi_l)
deallocate(config%jflde_l)
deallocate(config%mx_l)
deallocate(config%kxi_l)
deallocate(config%kxe_l)
deallocate(config%nm_l)
deallocate(config%ns_l)
deallocate(config%jsi_l)
deallocate(config%jse_l)
deallocate(config%kx_m)
deallocate(config%ky_m)
deallocate(config%jproc_m)

! close output file
close(unit=20)

end subroutine cleanup
