#' Display data.frame by using datatable function
#'
#' @param data displayed data frame.
#' @param scrollY Use page or scroller to display data. Default is TRUE.
#' @param pageLength the length of table. Default is 10 rows when scrollY is FALSE.
#' @export
show_db <- function(data, scrollY = TRUE, pageLength = 10){
  scrollY <- if(scrollY){
    scroller <- TRUE
    scrollY <- "40vh"
  }else{
    scroller <- FALSE
    scrollY <- NULL
  }
  DT::datatable(data, filter = 'top', extensions = c("Buttons", "ColReorder", "Scroller"),
                options = list(pageLength = pageLength, autoWidth = TRUE, rownames = FALSE,
                               dom = 'Bfrtip', buttons = I(c('colvis',"excel")), colReorder = TRUE,
                               deferRender = TRUE, scrollX = TRUE, scrollY = scrollY, scroller = scroller,
                               search = list(regex = TRUE, caseInsensitive = FALSE) )
  )
}
