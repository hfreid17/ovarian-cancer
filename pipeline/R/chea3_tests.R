{
 h

  df_test_forloop2 = df_test_forloop[-c(, 1)]

  list_dfs = list()

  for(i in names(df_test_forloop2)) {
    list_genes_column <- df_test_forloop2[[i]]
    char_array <- as.character(list_genes_column)
    result_inforloop = chea3::queryChea(geneset = char_array, set_name = i, n_results = "all", background = 20000)
    df_result1 = result_inforloop[[integrated]]
    df_result2 = df_result1[ , 10:11]
    list_dfs.append(df_result2)
  
  }
}


