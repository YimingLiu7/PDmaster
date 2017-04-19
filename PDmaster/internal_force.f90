!***************************************************************
!                    internal_force                            *
!          功能：                                              *
!***************************************************************    
subroutine internal_force()
    use scalars
    use serialArrays
    use globalParameters
    implicit none
    if (model=='PMB')then
       call model_bond_PMB()
    elseif(model=='LPS')then
       call model_state_LPS()
    elseif(model=='CDG')then
       call model_bond_CDG()
    else
       write(*,*)model,'材料模型类型出错！！！'
       stop
    end if
end subroutine internal_force
    
    
  



 