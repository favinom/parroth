[Mesh]
  [cyl2d_iga]
    type = FileMeshGenerator
    file = 'singleInter.xda'
  []
  [breakmesh]
    input = cyl2d_iga
    type = BreakMeshByBlockGenerator
    block_pairs = '2 4'
    split_interface = true
    #add_interface_on_two_sides = true
  []


[]





[Variables]
  [disp_x] []
  [disp_y] []
[]

[Kernels]
 [lex]
  type = LinearElasticity variable = disp_x disp_x = disp_x disp_y = disp_y component = 0
 []
 [ley]
  type = LinearElasticity variable = disp_y disp_x = disp_x disp_y = disp_y component = 1
 []
[]

[Materials]
 [./ElasticityMaterialProperties0] type = ElasticityMaterialProperties mu  = 1.0 lambda = 1.0 disp_x = disp_x disp_y = disp_y [../]
[]


[InterfaceKernels]
 [sigma_xx_jump_comp] type = ElasticityInterfaceKernel eps=1e-3 mu=0.001 lambda=0.01 variable = disp_x neighbor_var = disp_x
                 disp_x = disp_x disp_y = disp_y boundary = 'Block2_Block4' component = 0 jump = 1 []
 [sigma_yx_jump_comp] type = ElasticityInterfaceKernel eps=1e-3 mu=0.001 lambda=0.01 variable = disp_y neighbor_var = disp_y
                 disp_x = disp_x disp_y = disp_y boundary = 'Block2_Block4' component = 1 jump = 1 []

 [sigma_xx_mean_comp] type = ElasticityInterfaceKernel eps=1e-3 mu=0.001 lambda=0.01 variable = disp_x neighbor_var = disp_x disp_x = disp_x disp_y = disp_y boundary = 'Block2_Block4' component = 0 jump = 0 []
 [sigma_yx_mean_comp] type = ElasticityInterfaceKernel eps=1e-3 mu=0.001 lambda=0.01 variable = disp_y neighbor_var = disp_y disp_x = disp_x disp_y = disp_y boundary = 'Block2_Block4' component = 1 jump = 0 []
[]


[BCs]
 [bottom_x] type = DirichletBC value = 0 variable = disp_x boundary = bottom []
 [bottom_y] type = DirichletBC value = 0 variable = disp_y boundary = bottom []
 [top_x] type = DirichletBC value = 0.1 variable = disp_x boundary = top []
 [top_y] type = DirichletBC value = 0.1 variable = disp_y boundary = top []
[]

[Preconditioning]
  [smp]
    type = FDP  #SMP usando el jacobiano  - FDP usando una aproximacion del jacobiano
    full = true
  []
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  line_search = none
#  nl_rel_tol = 1e-10
#  nl_abs_tol = 1e-10

# petsc_options = '-ksp_view_pmat'
 petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
 petsc_options_value='  preonly   lu       NONZERO               mumps '

# [Quadrature] type = TRAP []

[]

[Outputs]
  file_base = inter_fracture_inclined
  print_linear_residuals = false
  exodus = true
#  xda = true
[]

