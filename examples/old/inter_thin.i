[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    nx = 20
    ny = 20
    dim = 2
  []
  [block1]
    type = SubdomainBoundingBoxGenerator
    block_id = 1
    bottom_left = '0.0 0.0 0.0'
    top_right =   '0.5 1.0 0.0'
    input = gen
  []

  [breakmesh]
    input = block1
    type = BreakMeshByBlockGenerator
    #block_pairs = '1 2'
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
 [./ElasticityMaterialProperties0] type = ElasticityMaterialProperties mu  = 1.0   lambda = 1.0 disp_x = disp_x disp_y = disp_y [../]
[]


[InterfaceKernels]
 [sigma_xx_jump] type = sigma_xx_jump eps=1e-3 mu=0.001 variable = disp_x neighbor_var = disp_x
                 disp_y = disp_y boundary = 'Block0_Block1' []
 [sigma_yx_jump] type = sigma_yx_jump eps=1e-3 mu=0.001 variable = disp_y neighbor_var = disp_y
                 boundary = 'Block0_Block1' []
 [sigma_xx_mean] type = sigma_xx_mean eps=1e-3 mu=0.001 variable = disp_x neighbor_var = disp_x boundary = 'Block0_Block1' []
 [sigma_yx_mean] type = sigma_yx_mean eps=1e-3 mu=0.001 variable = disp_y neighbor_var = disp_y disp_x = disp_x boundary = 'Block0_Block1' []
[]


[BCs]
 [left_x] type = DirichletBC value = 0 variable = disp_x boundary = left []
 [left_y] type = DirichletBC value = 0 variable = disp_y boundary = left []
 [right_x] type = DirichletBC value = 0.1 variable = disp_x boundary = right []
 [right_y] type = DirichletBC value = 0.1 variable = disp_y boundary = right []
[]

[Preconditioning]
  [smp]
    type = FDP
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
  file_base = solution2
  print_linear_residuals = false
  exodus = true
#  xda = true
[]

