library(shiny)
library(dplyr)
library(tidyr)

error_table <- readRDS("list_errors.RDS")
loglik_table <- readRDS("list_errors_loglik.RDS")
error_inla <- readRDS("list_inla_errors.RDS")
loglik_inla <- readRDS("list_errors_loglik_inla.RDS")

nu_val_loglik = c("0.1","0.4","0.6","0.8","1.2", "1.8", "2.2","2.8","3.1")


# GeomSpliViolin by jan-glx at https://stackoverflow.com/a/45614547

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

server <- function(input, output, session) {
  
  observe({
    mesh_size <- input$meshSize
    range_val <- input$rangeParameter
    norm_fun <- input$approxNorm
    trick <- input$useTrick
    m_order <- input$orderRat
    which_approx <- input$whichApprox
    includeINLA <- input$includeINLAcoverror
    
    color_plot_options <- c("steelblue", 
                    "limegreen", "red",
                    "purple")
    
    error_table_inla_cov <- as_tibble(error_inla[[mesh_size]][[range_val]])
    
    error_data_frame_inla <- error_table_inla_cov %>% select(nu, error_l2_inla) %>%
          rename("Covariance-based" = error_l2_inla) %>% mutate(Norm = "L2", Order = "0 (INLA)")
    tmp_tbl <- error_table_inla_cov %>% select(nu, error_sup_inla) %>%
          rename("Covariance-based" = error_sup_inla) %>% mutate(Norm = "Sup", Order = "0 (INLA)")
    
    error_data_frame_inla <- bind_rows(error_data_frame_inla, tmp_tbl)
    
    if(!is.null(m_order)){
      color_plot_used <- color_plot_options[as.numeric(m_order)]
      
      m_order_tmp <- as.character(m_order[1])
      error_table_cov <- as_tibble(error_table[[m_order_tmp]][[mesh_size]][[range_val]]) %>% 
         select(nu, error_l2_cov_trick, error_l2_op_trick) %>% mutate(Order = m_order_tmp) %>%
                   rename("Covariance-based" = error_l2_cov_trick, "Operator-based" = error_l2_op_trick) %>% mutate(Norm = "L2", Trick = "Yes")
      
      error_data_frame <- error_table_cov
      
      error_table_cov <- as_tibble(error_table[[m_order_tmp]][[mesh_size]][[range_val]]) %>% 
        select(nu, error_sup_cov_trick, error_sup_op_trick) %>% mutate(Order = m_order_tmp) %>% 
        rename("Covariance-based" = error_sup_cov_trick, "Operator-based" = error_sup_op_trick) %>% mutate(Norm = "Sup", Trick = "Yes")
      
      error_data_frame <- bind_rows(error_data_frame, error_table_cov)
      
      error_table_cov <- as_tibble(error_table[[m_order_tmp]][[mesh_size]][[range_val]]) %>% 
        select(nu, error_l2_cov_notrick, error_l2_op_notrick) %>% mutate(Order = m_order_tmp) %>% 
        rename("Covariance-based" = error_l2_cov_notrick, "Operator-based" = error_l2_op_notrick) %>% mutate(Norm = "L2", Trick = "No")
      
      error_data_frame <- bind_rows(error_data_frame, error_table_cov)
      
      error_table_cov <- as_tibble(error_table[[m_order_tmp]][[mesh_size]][[range_val]]) %>% 
        select(nu, error_sup_cov_notrick, error_sup_op_notrick) %>% mutate(Order = m_order_tmp) %>% 
        rename("Covariance-based" = error_sup_cov_notrick, "Operator-based" = error_sup_op_notrick) %>% mutate(Norm = "Sup", Trick = "No")
      
      error_data_frame <- bind_rows(error_data_frame, error_table_cov)
      
      if(length(m_order)>1){
        for(j in 2:length(m_order)){
          m_order_tmp <- as.character(m_order[j])
          tmp <- as_tibble(error_table[[m_order_tmp]][[mesh_size]][[range_val]]) %>% 
            select(nu, error_l2_cov_trick, error_l2_op_trick) %>% mutate(Order = m_order_tmp) %>% 
            rename("Covariance-based" = error_l2_cov_trick, "Operator-based" = error_l2_op_trick) %>% mutate(Norm = "L2", Trick = "Yes")
          error_data_frame <- bind_rows(error_data_frame, tmp)
          tmp <- as_tibble(error_table[[m_order_tmp]][[mesh_size]][[range_val]]) %>% 
            select(nu, error_sup_cov_trick, error_sup_op_trick) %>% mutate(Order = m_order_tmp) %>% 
            rename("Covariance-based" = error_sup_cov_trick, "Operator-based" = error_sup_op_trick) %>% mutate(Norm = "Sup", Trick = "Yes")
          error_data_frame <- bind_rows(error_data_frame, tmp)
          tmp <- as_tibble(error_table[[m_order_tmp]][[mesh_size]][[range_val]]) %>% 
            select(nu, error_l2_cov_notrick, error_l2_op_notrick) %>% mutate(Order = m_order_tmp) %>% 
            rename("Covariance-based" = error_l2_cov_notrick, "Operator-based" = error_l2_op_notrick) %>% mutate(Norm = "L2", Trick = "No")
          error_data_frame <- bind_rows(error_data_frame, tmp)
          tmp <- as_tibble(error_table[[m_order_tmp]][[mesh_size]][[range_val]]) %>% 
            select(nu, error_sup_cov_notrick, error_sup_op_notrick) %>% mutate(Order = m_order_tmp) %>% 
            rename("Covariance-based" = error_sup_cov_notrick, "Operator-based" = error_sup_op_notrick) %>% mutate(Norm = "Sup", Trick = "No")
          error_data_frame <- bind_rows(error_data_frame, tmp)
        }
      }
      
      error_data_frame <- error_data_frame %>% filter(Norm == norm_fun)
      error_data_frame_inla <- error_data_frame_inla %>% filter(Norm == norm_fun)
      
      error_data_frame <- error_data_frame %>% tidyr::pivot_longer(cols = ends_with("-based"), names_to = "Type", values_to = "error")
      error_data_frame_inla <- error_data_frame_inla %>% tidyr::pivot_longer(cols = ends_with("-based"), names_to = "Type", values_to = "error")
      
      if(trick){
        error_data_frame <- error_data_frame %>% filter(Trick == "Yes")
      } else{
        error_data_frame <- error_data_frame %>% filter(Trick == "No")
      }
      
      if(!is.null(input$whichApprox)){
        error_data_frame <- error_data_frame %>% filter(Type %in% which_approx)
      } else{
        updateCheckboxGroupInput(session, "whichApprox", selected= c("Covariance-based", "Operator-based"))
      }

      if(includeINLA){
        color_plot_used <- c("black", color_plot_used)
      }
      
      
      fig <- ggplot(error_data_frame, aes(x = nu, y = error, linetype=Type, color = Order,
                                   group = interaction(Type, Order))) + geom_line() +
        scale_color_manual(values = color_plot_used)
      
      if(includeINLA){
        fig <- fig + geom_line(data = error_data_frame_inla, aes(y = error)) 
      }
      
      y_label_cov <- ifelse(norm_fun == "L2", "Error in L2-norm", "Error in Sup-norm")
      
      fig <- fig +labs(y= y_label_cov, x = "\u028b (smoothness parameter)")
      
      if(input$logScaleCoverror){
        fig <- fig + scale_y_log10()
      }
      
      output$downloadPlotCov <- downloadHandler(
        filename = function() { paste0("Cov_error_plot.", input$fileExtensionCov) },
        content = function(file) {
          ggsave(file, plot = fig, device = input$fileExtensionCov)
        }
      )
      
      output$downloadDataFrameCov <- downloadHandler(
        filename = function() { paste0("Cov_error_data_frame_range_",range_val,".RDS") },
        content = function(file) {
          saveRDS(bind_rows(error_data_frame,error_data_frame_inla), file = file)
        }
      )
      
      fig_plotly <- ggplotly(fig)
      
      
    } else{
      m_order = c(2)
      updateCheckboxGroupInput(session, "orderRat", selected= c("2"))
    }
    
    output$raterrors <- renderPlotly({
      
      fig_plotly
      
    })
    
    
  })
  
  
  observe({
    range_val <- input$rangeParameterLoglik
    m_order <- input$orderRatLoglik
    sigma.e <- input$measError
    
    error_table_loglik_cov <- as_tibble(loglik_table[[m_order]][[range_val]][[sigma.e]][["covariance_based"]])
    error_table_loglik_op <- as_tibble(loglik_table[[m_order]][[range_val]][[sigma.e]][["operator_based"]])
    
    error_table_inla <- as_tibble(loglik_inla[[range_val]][[sigma.e]][["inla"]]) 
    
    if(input$useAbsRelError){
      error_table_loglik_cov <- abs(error_table_loglik_cov)
      error_table_loglik_op <- abs(error_table_loglik_op)
      error_table_inla <- abs(error_table_inla)
    }
    
    error_table_loglik_cov <- error_table_loglik_cov %>% mutate(Type = "Covariance-based") %>% select(-Rel.error)
    error_table_loglik_op <- error_table_loglik_op %>% mutate(Type = "Operator-based") %>% select(-Rel.error)
    error_table_inla <- error_table_inla %>% mutate(Type = "INLA") %>% select(-Rel.error)

    lik_table <- bind_rows(error_table_loglik_cov, error_table_loglik_op)
    
    lik_table <- lik_table %>% tidyr::pivot_longer(
      cols = starts_with("nu =  "), 
      names_to = "nu", 
      values_to = "lik_error", 
      names_prefix = "nu =  ")
    
    error_table_inla <- error_table_inla  %>% tidyr::pivot_longer(
      cols = starts_with("nu =  "), 
      names_to = "nu", 
      values_to = "lik_error", 
      names_prefix = "nu =  ")
    
    nu_loglik <- input$smoothloglik
    
    if(!is.null(nu_loglik)){
    
      lik_table <- lik_table %>% filter(nu %in% nu_loglik)
      error_table_inla <- error_table_inla %>% filter(nu %in% nu_loglik)
      
      if(input$whichComparisonsLoglik == "Operator-based"){
        lik_table <- lik_table %>% filter(Type == "Operator-based")
      } else if(input$whichComparisonsLoglik == "Covariance-based"){
        lik_table <- lik_table %>% filter(Type == "Covariance-based")
      }
      
    fig <- ggplot(lik_table, aes(x = nu, y = lik_error, fill = Type))
    
    if(input$includeINLA){
      fig <- fig +geom_violin(data = error_table_inla, aes(y = lik_error), scale="width") 
    }
    
    if(input$whichComparisonsLoglik == "Both"){
      fig <- fig + geom_split_violin(scale="width")
    } else{
      fig <- fig +geom_violin(data = lik_table, aes(y = lik_error), scale="width") 
    }
    
    
    
    fig <- fig +labs(y= "Loglikelihood relative error",  x = expression(nu~"(smoothness parameter)"))
    
    
    output$downloadPlotLogLik <- downloadHandler(
      filename = function() { paste0("LogLik_error_plot.", input$fileExtensionLogLik) },
      content = function(file) {
        ggsave(file, plot = fig, device = input$fileExtensionLogLik)
      }
    )
    
    output$downloadDataFrameViolin <- downloadHandler(
      filename = function() { paste0("Violin_loglik_data_frame_range_",range_val,".RDS") },
      content = function(file) {
        saveRDS(bind_rows(lik_table,error_table_inla), file = file)
      }
    )
    
    } else{
      updateCheckboxGroupInput(session, "smoothloglik", selected= nu_val_loglik)
      nu_loglik <- nu_val_loglik
    }

    output$loglikerrorsfig <- renderPlot({

      fig

    })
  })
  
  
  observe({
    range_val <- input$rangeParameterLoglikMean
    m_order <- input$orderRatLoglikMean
    sigma.e <- input$measErrorMean
    nu_loglik <- input$smoothloglikMean
    
    if(!is.null(m_order)){
      m_order <- as.numeric(m_order)
      m_order_tmp <- m_order[1]
      error_table_loglik_cov <- as_tibble(loglik_table[[m_order_tmp]][[range_val]][[sigma.e]][["covariance_based"]]) %>% mutate(Order = m_order_tmp)
      error_table_loglik_op <- as_tibble(loglik_table[[m_order_tmp]][[range_val]][[sigma.e]][["operator_based"]]) %>% mutate(Order = m_order_tmp)
      if(length(m_order)>1){
        for(j in 2:length(m_order)){
          m_order_tmp <- m_order[j]
          tmp <- as_tibble(loglik_table[[m_order_tmp]][[range_val]][[sigma.e]][["covariance_based"]]) %>% mutate(Order = m_order_tmp)
          error_table_loglik_cov <- bind_rows(error_table_loglik_cov, tmp)
          tmp <- as_tibble(loglik_table[[m_order_tmp]][[range_val]][[sigma.e]][["operator_based"]]) %>% mutate(Order = m_order_tmp)
          error_table_loglik_op <- bind_rows(error_table_loglik_op, tmp)
        }
      }
      
      
      
      
      
      error_table_inla <- as_tibble(loglik_inla[[range_val]][[sigma.e]][["inla"]]) 
      
      if(input$useAbsRelErrorMean){
        error_table_loglik_cov <- abs(error_table_loglik_cov)
        error_table_loglik_op <- abs(error_table_loglik_op)
        error_table_inla <- abs(error_table_inla)
      }
      
      error_table_loglik_cov <- error_table_loglik_cov %>% select(-Rel.error)
      error_table_loglik_op <- error_table_loglik_op %>% select(-Rel.error) 
      error_table_inla <- error_table_inla %>% select(-Rel.error) %>% summarise(across(everything(), mean), .groups = 'drop')
      
      error_table_loglik_cov <- error_table_loglik_cov %>% mutate(Type = "Covariance-based") 
      error_table_loglik_op <- error_table_loglik_op %>% mutate(Type = "Operator-based") 
      error_table_inla <- error_table_inla %>% mutate(Type = "INLA")
      
      lik_table <- bind_rows(error_table_loglik_cov, error_table_loglik_op) %>% group_by(Type, Order) %>% summarise(across(everything(), mean), .groups = 'drop')
      
      lik_table <- lik_table %>% tidyr::pivot_longer(
        cols = starts_with("nu =  "), 
        names_to = "nu", 
        values_to = "lik_error", 
        names_prefix = "nu =  ")
      
      error_table_inla <- error_table_inla  %>% tidyr::pivot_longer(
        cols = starts_with("nu =  "), 
        names_to = "nu", 
        values_to = "lik_error", 
        names_prefix = "nu =  ")
      
      if(!is.null(nu_loglik)){
        
        lik_table <- lik_table %>% filter(nu %in% nu_loglik) %>% mutate(nu = as.numeric(nu), Order = as.character(Order)) 
        error_table_inla <- error_table_inla %>% filter(nu %in% nu_loglik) %>% mutate(nu = as.numeric(nu), Order = as.character(0))
        
        if(input$whichComparisonsLoglikMean == "Operator-based"){
          lik_table <- lik_table %>% filter(Type == "Operator-based")
        } else if(input$whichComparisonsLoglikMean == "Covariance-based"){
          lik_table <- lik_table %>% filter(Type == "Covariance-based")
        }
        
        fig <- ggplot(lik_table, aes(x = nu, y = lik_error, linetype=Type, color = Order,
                                     group = interaction(Type, Order))) + geom_line()
        
        if(input$includeINLAmean){
          fig <- fig +geom_line(data = error_table_inla, aes(y = lik_error)) 
        }
        
        fig <- fig +labs(y= "Loglikelihood relative error (mean)", x = "\u028b (smoothness parameter)") 
        
        if(input$logScaleMean){
          fig <- fig + scale_y_log10()
        }
        
        fig_plotly <- ggplotly(fig)
        
        output$downloadPlotLogLikMean <- downloadHandler(
          filename = function() { paste0("LogLikMean_error_plot.", input$fileExtensionLogLikMean) },
          content = function(file) {
            ggsave(file, plot = fig, device = input$fileExtensionLogLikMean)
          }
        )
        
        output$downloadDataFrameMean <- downloadHandler(
          filename = function() { paste0("Mean_loglik_data_frame_range_",range_val,".RDS") },
          content = function(file) {
            saveRDS(bind_rows(lik_table,error_table_inla), file = file)
          }
        )
        
      } else{
        updateCheckboxGroupInput(session, "smoothloglikMean", selected= nu_val_loglik)
        nu_loglik <- nu_val_loglik
      }
      
      output$loglikerrorsfigMean <- renderPlotly({
        
        fig_plotly
        
      })
      

    } else{
      m_order = c(2)
      updateCheckboxGroupInput(session, "orderRatLoglikMean", selected= c(2))
    }

    
  })
  
}