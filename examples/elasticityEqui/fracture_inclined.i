[Mesh]
  [cyl2d_iga]
    type = FileMeshGenerator
    file = 'singleEqui.xda'
  []
[]

[Variables]
  [disp_x] []
  [disp_y] []
[]

[Kernels]
 [lex]
  type = LinearElasticity variable = disp_x disp_x = disp_x disp_y = disp_y  component = 0
 []
 [ley]
  type = LinearElasticity variable = disp_y disp_x = disp_x disp_y = disp_y  component = 1
 []

[]

[Materials]
 [./ElasticityMaterialProperties0] type = ElasticityMaterialProperties mu  = 1.0   lambda = 1.0  disp_x = disp_x disp_y = disp_y  block = '1 2 4' [../]
 [./ElasticityMaterialProperties1] type = ElasticityMaterialProperties mu  = 0.001 lambda = 0.01 disp_x = disp_x disp_y = disp_y  block = 3       [../]
[]


[BCs]
 [bottom_x] type = DirichletBC value = 0 variable = disp_x boundary = bottom []
 [bottom_y] type = DirichletBC value = 0 variable = disp_y boundary = bottom []
 [top_x] type = DirichletBC value = 0.1 variable = disp_x boundary = top []
 [top_y] type = DirichletBC value = 0.1 variable = disp_y boundary = top []
[]


[Preconditioning]
  [smp]
    type = SMP
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
  file_base = equi_fracture_inclined
  print_linear_residuals = false
  exodus = true
#  xda = true
[]

