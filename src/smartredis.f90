
MODULE mod_smartredis
  use mod_precision, only: rp
  USE smartredis_client, ONLY: CLIENT_TYPE
  IMPLICIT NONE
  TYPE(CLIENT_TYPE)   :: client

  CONTAINS

  SUBROUTINE InitSmartRedis()
    IMPLICIT NONE
    integer :: error

    error = client%Initialize(.false.)

  END SUBROUTINE InitSmartRedis


  SUBROUTINE ExchangeDataSmartRedis(vel1, vel2, tauw1, tauw2)
    IMPLICIT NONE
    real(rp), intent(in), dimension(0:,0:) :: vel1, vel2
    real(rp), intent(out), dimension(0:,0:) :: tauw1, tauw2

    logical :: exists
    integer :: error

    error = client%put_tensor("train000_state", vel1, shape(vel1))

    error = client%put_tensor("train000_reward", (/0.1/), (/1/))

    exists = .false.
    print *, "polling for train000_actions"
    error = client%poll_tensor("train000_actions", 10, 100, exists) ! 10ms*1000 = 10s
    print *, "polling for train000_actions done"

    if (.not.exists) then
      print*, "poll error"
    else
      error = client%unpack_tensor("train000_actions", tauw1, shape(tauw1))
      error = Client%delete_tensor("train000_actions")
    end if

  END SUBROUTINE ExchangeDataSmartRedis


  SUBROUTINE FinalizeSmartRedis()
    IMPLICIT NONE
    integer :: error
    error = client%destructor()
  END SUBROUTINE FinalizeSmartRedis

END MODULE mod_smartredis
