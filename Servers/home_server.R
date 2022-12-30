# Channge the tab when a button is pressed. 
observeEvent(input$tut_switch, {
  updateTabsetPanel(inputId = "main_tab", selected = "Tutorial")
})
observeEvent(input$data_switch, {
  updateTabsetPanel(inputId = "main_tab", selected = "Data Upload")
})
observeEvent(input$ex_switch, {
  updateRadioButtons(inputId = "data_type", selected = "ed")
  
  updateTabsetPanel(inputId = "main_tab", selected = "Visualization")
})