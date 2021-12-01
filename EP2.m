0;
#EP2
#Nome : Bruno Hideki Akamine
#Nusp : 11796322

function [ind x d] = simplex_res(A,b,c,m,n)
  m_backup = m;
  inv = [eye(m)];
  for i = 1 : m
    if b(i) < 0
      b(i) = - b(i);
      A(i, :) = - A(i, :);
    endif
  endfor
  A = [A inv];    #Note que a inversa da matriz identidade é a própria identidade
  colunas_basicas = (n + 1) : (m + n);
  colunas_Nbasicas = 1 : n;
  c_Aux = [zeros(n, 1)', ones(m, 1)'];
  x = [zeros(n, 1)', b];
  printf("----------      Fase 1      ----------\n")
  iter = 0;
  while true
    c_Aux_b = c_Aux(colunas_basicas);
    printf("***********Iteração %d***************\n", iter)
    printf("Indices das variáveis básicas:(indice - valor associado)\n")
    for i = 1:length(colunas_basicas)
      printf("%d  %.05f\n", colunas_basicas(i), x(colunas_basicas(i)))
    endfor
    printf("\nValor da função objetiva: %.05f\n\n", dot(x, c_Aux))
    printf("Matriz inversa:\n")
    disp(inv)
    printf("\n")
    p = c_Aux_b * inv;
    p = p';
    printf("Vetor p:\n")
    disp(p)
    enter = m + n + 1;
    printf("Custos reduzidos: \n")
    for i = 1:length(colunas_Nbasicas)
      cj = c_Aux(colunas_Nbasicas(i)) - dot(p' , A(:, colunas_Nbasicas(i)));
      printf("%d  %.05f\n", colunas_Nbasicas(i), cj)
      if cj < 0 && enter > colunas_Nbasicas(i) 
        enter = colunas_Nbasicas(i);
      endif
    endfor
    if enter < m + n + 1
      printf("\n Entra na base: %d\n\n", enter)
    else
      break;
    endif
    u = inv * A(:, enter);
    exit = n + m + 1;
    printf("Direção: \n")
    for i = 1:length(u)
      printf("%d  %.05f\n", colunas_basicas(i), u(i))
    endfor 
    teta = inf;
    for i = 1:length(u)
        if u(i) > 0
          prov = x(colunas_basicas(i))/u(i);
          if teta == prov
            if exit > colunas_basicas(i);
              exit = colunas_basicas(i);
            endif
          endif
          if teta > prov
            teta = prov;
            exit = colunas_basicas(i);
          endif
        endif
    endfor
    printf("\nTheta*\n  %.05f\n", teta)
    printf("\nSai da base : %d\n\n", exit)
    for i = 1:length(colunas_basicas)
        if colunas_basicas(i) == exit
          index = i;
          x(colunas_basicas(i)) = 0;
          colunas_basicas(i) = enter;
        else
          x(colunas_basicas(i)) = x(colunas_basicas(i)) - teta * u(i);
        endif
    endfor
    for i = 1 : length(colunas_Nbasicas)
        if colunas_Nbasicas(i) == enter
          x(colunas_Nbasicas(i)) = teta;
          colunas_Nbasicas(i) = exit;
        endif
    endfor
    iter = iter + 1;
    exitrow = inv(index, :);
    for i = 1:m
      if i != index
        inv(i,:) = inv(i, :) + (exitrow * (-u(i)/u(index)));
       else
        inv(i, :) = inv(i, :) / u(index);
      endif
    endfor
  endwhile
  if dot(x, c_Aux)> 0
    ind = 1;
    d = NaN;
    x = NaN;
    printf("Problema inviável.\n")
    return;
  endif
  for i = 1 : length(colunas_basicas)
    if i > length(colunas_basicas)
      break;
    endif
    if colunas_basicas(i) > n
       printf("\n    -Retirando a variável artificial %d da base-\n", colunas_basicas(i))
       enter = n + 1;
       exit = colunas_basicas(i);
       prov = inv * A(:, 1:n);
       printf("\nLinha correspondente a variavel %d da matriz inv * A:\n", colunas_basicas(i))
       for j = 1 : n
           if prov(i, j) != 0 && enter > n
               enter = j;
           endif
           printf("%.05f ", prov(i, j))
       endfor
       printf("\n")
       if enter < n + 1
          printf("\nVariável %d entra na base.\n", enter)
          u = inv * A(:, enter);
          teta = x(exit)/u(i);
          for k = 1:length(colunas_basicas)
            if colunas_basicas(k) == exit
              index = k;
              x(colunas_basicas(k)) = 0;
              colunas_basicas(k) = enter;
            else
              x(colunas_basicas(k)) = x(colunas_basicas(k)) - teta * u(k);
            endif
          endfor
          for k = 1 : length(colunas_Nbasicas)
              if colunas_Nbasicas(k) == enter
                x(colunas_Nbasicas(k)) = teta;
                colunas_Nbasicas(k) = exit;
              endif
          endfor
          exitrow = inv(index, :);
          for k = 1:m
            if k != index
              inv(k, :) = inv(k, :) + (exitrow * (-u(k)/u(index)));
            else
              inv(k, :) = inv(k, :) / u(index);
            endif
          endfor
       else
          printf("\nLinha %d da matriz A corresponde a uma restrição redundante, portanto vamos retirá-la.\n", i)
          printf("Retiramos também a variável %d das variáveis básicas.\n", colunas_basicas(i))
           A = [A(1: i - 1, :); A(i + 1: m , :)];
           b = [b(1:i - 1) b(i + 1: m)];
           inv = [inv(:, 1:i - 1) inv(:, i + 1: m)];
           inv = [inv(1 : i - 1, :);inv(i + 1 : m, :)];
           colunas_Nbasicas = [colunas_Nbasicas colunas_basicas(i)];
           colunas_basicas = [colunas_basicas(1:i - 1) colunas_basicas(i + 1:m)];
           m = m - 1;
       endif
       printf("Matriz inversa:\n")
       disp(inv)
       printf("\n")
    endif
  endfor
  A = A(:, 1 : n);
  colunas_Nbasicas = setdiff(colunas_Nbasicas, n + 1 : n + m_backup);
  x = x(1 : n);
  printf("----------      Fase 2      ----------\n")
  iter = 0;
  while true
    cb = c(colunas_basicas);
    printf("***********Iteração %d***************\n", iter)
    printf("Indices das variáveis básicas:(indice - valor associado)\n")
    for i = 1:length(colunas_basicas)
      printf("%d  %.05f\n", colunas_basicas(i), x(colunas_basicas(i)))
    endfor
    printf("\nValor da função objetiva: %.05f\n\n", dot(x, c))
    printf("Matriz inversa: \n")
    disp(inv)
    printf("\n")
    p = cb * inv;
    p = p';
    printf("Vetor p:\n")
    disp(p)
    printf("\n")
    enter = n + 1;
    printf("Custos reduzidos: \n")
    for i = 1:length(colunas_Nbasicas)
      cj = c(colunas_Nbasicas(i)) - dot(p' , A(:, colunas_Nbasicas(i)));
      printf("%d  %.05f\n", colunas_Nbasicas(i), cj)
      if cj < 0 && enter > colunas_Nbasicas(i) 
        enter = colunas_Nbasicas(i);
      endif
    endfor
    if enter < n + 1
      printf("\n Entra na base: %d\n\n", enter)
    else
      ind = 0;
      printf("\nSolução ótima encontrada com custo %.05f:\n", dot(x, c))
      for i = 1: length(x)
        printf("%d  %.05f\n", i, x(i))
      endfor
      d = NaN;
      break;
    endif
    u = inv * A(:, enter);
    exit = n + 1;
    printf("Direção: \n")
    for i = 1:length(u)
      printf("%d  %.05f\n", colunas_basicas(i), u(i))
    endfor 
    for i = 1:length(u)
      if u(i) > 0
        exit = colunas_basicas(i);
        break;
      endif
    endfor
    if exit < n + 1
      teta = Inf;
      for i = 1:length(u)
        if u(i) > 0
          prov = x(colunas_basicas(i))/u(i);
          if teta == prov
            if exit > colunas_basicas(i);
              exit = colunas_basicas(i);
            endif
          endif
          if teta > prov
            teta = prov;
            exit = colunas_basicas(i);
          endif
        endif
      endfor
      printf("\nTheta*\n  %.05f\n", teta)
      printf("\nSai da base : %d\n\n", exit)
      for i = 1:length(colunas_basicas)
        if colunas_basicas(i) == exit
          index = i;
          x(colunas_basicas(i)) = 0;
          colunas_basicas(i) = enter;
        else
          x(colunas_basicas(i)) = x(colunas_basicas(i)) - teta * u(i);
        endif
      endfor
      for i = 1 : length(colunas_Nbasicas)
        if colunas_Nbasicas(i) == enter
          x(colunas_Nbasicas(i)) = teta;
          colunas_Nbasicas(i) = exit;
        endif
      endfor
    else
      printf("Problema ilimitado\n")
      d = u;
      ind = -1;
      x = NaN;
      break;
    endif
    iter = iter + 1;
    exitrow = inv(index, :);
    for i = 1:m
      if i != index
        inv(i,:) = inv(i, :) + (exitrow * (-u(i)/u(index)));
       else
        inv(i, :) = inv(i, :) / u(index);
      endif
    endfor
  endwhile
  
endfunction

function [ind x d] = simplex_tab(tableau)
  c = tableau(1, :);
  c = c(2:length(c));
  b = tableau(:, 1);
  b = b(2:length(b));
  A = tableau(2:size(tableau, 1) , 2: size(tableau, 2));
  m = size(tableau, 1) - 1;
  m_backup = m;
  n = size(tableau, 2) - 1;
  linha_zero = [-sum(b), -ones(m ,1)' * A, zeros(m, 1)'];
  artificiais = eye(m);
  colunas_Nbasicas = 1:n;
  colunas_basicas = n + 1: n + m;
  tableau = tableau(2:size(tableau, 1), :);
  tableau = [tableau artificiais];
  tableau = [linha_zero ;tableau];
  iter = 0;
  printf("----------      Fase 1      ----------\n")
  while true
    printf("***********Iteração %d***************\n", iter)
    printf("Tableau:\n")
    disp(tableau)
    printf("\n")
    enter = n + m + 1;
    for i = colunas_Nbasicas
      if tableau(1, i + 1) < 0 && enter > i
        enter = i;
      endif
    endfor
    if enter < n + m + 1
      printf("\n Entra na base : %d\n\n", enter)
    else
      break;
    endif
    exit = n + m + 1;
    teta = Inf;
    for i = 1 : length(colunas_basicas)
      if tableau(i + 1, enter + 1) > 0
        prov = tableau(i + 1, 1)/tableau(i + 1, enter + 1);
        if teta == prov && exit > colunas_basicas(i)
          exit = i;
        endif
        if teta > prov
          teta = prov;
          exit = i;
        endif
      endif
    endfor
    printf("\nTheta*\n%.05f\n", teta)
    printf("\nSai da base : %d\n\n", colunas_basicas(exit))
    colunas_basicas(exit) = enter;
    iter = iter + 1;
    exitrow = tableau(exit + 1, :);
    u = tableau(:, colunas_basicas(exit) + 1);
    for i = 1:size(tableau, 1)
      if i != exit + 1
        tableau(i, :) = tableau(i, :) + ((-u(i)/u(exit + 1)) .* exitrow);
       else
        tableau(i, :) = tableau(i, :) ./ tableau(i , enter + 1);
      endif
    endfor
    colunas_Nbasicas = setdiff(1:size(tableau, 2) -1,colunas_basicas);
  endwhile
  if -tableau(1,1) > 0
    ind = 1;
    printf("Problema inviável.\n")
    return;
  endif
  for i = 1:length(colunas_basicas)
    if i > length(colunas_basicas)
      break;
    endif
    if colunas_basicas(i) > n
      printf("\n    -Retirando a variável artificial %d da base-\n", colunas_basicas(i))
      printf("\nTableau:\n")
      disp(tableau)
      printf("\n")
      enter = n + 1;
      exit = i;
      for j = 2 : n + 1
        if tableau(i + 1, j) != 0
          enter = j;
        endif
      endfor
      if enter > n
        printf("\nLinha %d do tableau corresponde a uma restrição redundante, portanto vamos retirá-la.\n", i + 1)
        printf("Retiramos também a variável %d das variáveis básicas.\n", colunas_basicas(i))
        colunas_Nbasicas = [colunas_Nbasicas colunas_basicas(i)];
        colunas_basicas = [colunas_basicas(1:i - 1) colunas_basicas(i + 1:m)];
        tableau = [tableau(1 : i , :) ;tableau(i + 2 : m + 1, :)];
        m = m - 1;
      else
        printf("\nVariável %d entra na base.\n", enter - 1)
        exitrow = tableau(exit + 1, :);
        u = tableau(:, colunas_basicas(exit) + 1);
        for i = 1:size(tableau, 1)
          if i != exit + 1
            tableau(i, :) = tableau(i, :) + ((-u(i)/u(exit + 1)) .* exitrow);
          else
            tableau(i, :) = tableau(i, :) ./ tableau(i , enter);
          endif
        endfor
        colunas_basicas(exit) = enter - 1;
        colunas_Nbasicas = setdiff(colunas_Nbasicas, n + 1 : n + m_backup);
      endif
    endif
  endfor
  printf("\nTableau final da fase 1:\n")
  disp(tableau)
  printf("\n")
  tableau = tableau(:, 1: n + 1);
  iter = 0;
  colunas_Nbasicas = setdiff(colunas_Nbasicas, n + 1 : n + m_backup);
  x = zeros(n ,1)';
  for i = 1 : length(colunas_basicas)
    x(colunas_basicas(i)) = tableau(i + 1, 1);
  endfor
  cb = c(colunas_basicas);
  linha_zero = [ -dot(x, c) c - cb * tableau(2 : m + 1, 2 : n + 1)];
  tableau = tableau(2 : m + 1, :);
  tableau = [linha_zero ; tableau];
  printf("----------      Fase 2      ----------\n")
  while true
    colunas_Nbasicas = setdiff(1:size(tableau, 2) -1,colunas_basicas);
    printf("***********Iteração %d***************\n", iter)
    printf("\nTableau:\n")
    disp(tableau)
    enter = n + 1;
    for i = 1 : n
      if tableau(1, i + 1) < 0
        enter = i;
        break;
      endif
    endfor
     if enter <= n
      printf("\n Entra na base: %d\n\n", enter)
    else
      ind = 0;
      x = [];
      printf("\nSolução ótima encontrada com custo %.05f:\n", -tableau(1, 1))
      for i = 1: size(tableau, 2) - 1
        x = zeros(size(tableau,2) - 1, 1);
        x = x';
        d = NaN;
        for i = 1:length(colunas_basicas)
          x(colunas_basicas(i)) = tableau(i + 1,1);
        endfor
      endfor
      for i = 1:length(x)
        printf("%d %.05f\n", i, x(i))
      endfor
      break;
    endif
     exit = n + 1;
    for i = 1:length(colunas_basicas)
      if tableau(i + 1, enter + 1) > 0
        exit = i;
        break;
      endif
    endfor
    if exit <= n
      teta = Inf;
      for i = 1:length(colunas_basicas)
        if tableau(i + 1, enter + 1) > 0
          prov = tableau(i + 1, 1)/tableau(i + 1, enter + 1);
          if teta > prov
            teta = prov;
            exit = i;
          endif
        endif
      endfor
      printf("\nTheta*\n%.05f\n", teta)
      printf("\nSai da base : %d\n\n", colunas_basicas(exit))
      colunas_basicas(exit) = enter;
    else
      printf("Problema ilimitado\n")
      d = tableau(2 : m + 1, enter + 1);
      x = NaN;
      ind = -1;
      break;
    endif
    iter = iter + 1;
    exitrow = tableau(exit + 1, :);
    u = tableau(:, colunas_basicas(exit) + 1);
    for i = 1:size(tableau, 1)
      if i != exit + 1
        tableau(i, :) = tableau(i, :) + ((-u(i)/u(exit + 1)) .* exitrow);
       else
        tableau(i, :) = tableau(i, :) ./ tableau(i , enter + 1);
      endif
    endfor
  endwhile
endfunction
#{
A = [1 3 0 4 1; 1 2 0 -3 1; -1 -4 3 0 0];
b = [2 2 1];
m = 3;
n = 5;
c = [2 3 3 1 -2];
tableau = [b' A];
tableau = [0 c; tableau]
simplex_res(A,b,c,m,n)
#simplex_tab(tableau)
#}

#{
A = [1 2 3 0; -1 2 6 0; 0 4 9 0; 0 0 3 1];
b = [3 2 5 1];
m = 4;
n = 4;
c = [1 1 1 0];
tableau = [b' A];
tableau = [0 c; tableau]
#[ind x d] = simplex_res(A,b,c,m,n)
[ind x d] = simplex_tab(tableau)
#}
#{
m = 2;
n = 4;
c = [-2 -1 0 0];
b = [1 2];
A = [-1 1 1 0; 1 -2 0 1];
#[ind x d] = simplex_res(A,b,c,m,n)
tableau = [b' A];
tableau = [0 c; tableau]
[ind x d] = simplex_tab(tableau)
#}
#{
m = 3; 
n = 4;
A = [1 1 1 1 ; 7 5 3 2 ; 3 5 10 15 ];
b = [15, 120, 100];
c = [-4 -5 -9 -11];
tableau = [b' A];
tableau = [0 c; tableau]
#[ind x d] = simplex_res(A,b,c,m,n)
simplex_tab(tableau)
#}
#{
A = [1 1 1 1 5; 1 1 1 1 8; 1 1 1 1 6];
c = [1 1 1 1 1];
b = [6 9 7];
m = 3;
n = 5;
[ind x d] = simplex_res(A,b,c,m,n)
#}