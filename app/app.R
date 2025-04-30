#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# 
# Copyright © 2025, Umeå Plant Science Center, Swedish University Of Agricultural Sciences, Umeå, Sweden
# All rights reserved.
# 
# Developers: Adrien Heymans
# 
# Redistribution and use in source and binary forms, with or without modification, are permitted under the GNU General Public License v3 and provided that the following conditions are met:
#   
#   1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# 
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# 
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
# 
# Disclaimer
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# You should have received the GNU GENERAL PUBLIC LICENSE v3 with this file in license.txt but can also be found at http://www.gnu.org/licenses/gpl-3.0.en.html
# 
# NOTE: The GPL.v3 license requires that all derivative work is distributed under the same license. That means that if you use this source code in any other program, you can only distribute that program with the full source code included and licensed under a GPL license.



library(shiny)
library(tidyverse)
source("./R/prep_geo.R")
source("./R/get_root_section.R")
source("./R/write_geo.R")
`%!in%` <- compose(`!`, `%in%`)

# Anatomy to mesh
ui <- fluidPage(
  titlePanel("Anatomy Deck"),
  
  sidebarLayout(
    sidebarPanel(
        h3('From ".xml" file using CellSeT format'),
        
        radioButtons(
          "input_mode", "Select source:",
          choices = c("Upload XML", "Upload ROIs", "Built-in anatomies"),
          selected = "Upload XML"
        ),
        
        conditionalPanel(
          condition = "input.input_mode == 'Upload XML'",
          fileInput("xml_file", "Upload Root XML File", accept = ".xml"),
          numericInput(inputId = "dpi", label = "Enter image scale (micron per pixel)", 
                       value = 2)
        ),
        
        conditionalPanel(
          condition = "input.input_mode == 'Upload ROIs'",
          fileInput("all_roi_files", "Upload ROI Files", accept = ".roi", multiple=T),
          numericInput(inputId = "dpi", label = "Enter image scale (micron per pixel)", 
                       value = 2)
        ),
        
        conditionalPanel(
          condition = "input.input_mode == 'Built-in anatomies'",
          selectInput("choice", "Choose an anatomical sample:",
                      choices = c("root cross section","root tip", "apical meristem", "embryon", "pavement cells")                      )
        ),
        
        
        sliderInput(inputId="thickness", label= "Cell wall thickness", 0, 2,
                    0.3, 0.1),
        sliderInput(inputId="corner", label= "Cell corner smoothing", 0.1, 5,
                    0.5, 0.1),
        h5("Once the original anatomy is loaded:"),
        actionButton("prep_btn", "Prep GEO"),
        h5("Save with GEO extension if you intend to mesh it with GSMH"),
      
      # Save
      radioButtons("export_type", "Export as:", choices = c("GEO", "PNG", "SVG", "CSV")),
      downloadButton("download_btn", "Download"),
      # 
      

    ),
    
    mainPanel(
      h4("Original Anatomy"),
      plotOutput("anatomy_plot"),
      
      h4("Preprocessed Geometry"),
      plotOutput("geo_plot")
    )
  )
)

# Define server logic required to make geo file
server <- function(input, output, session) {
  
  # Read XML into root section

  root_data <- reactive({
    req(input$dpi)
    if (input$input_mode == "Upload XML") {
      req(input$xml_file)
      section = get_root_section(input$xml_file$datapath)%>%
        mutate(x = x/input$dpi, y = -y/input$dpi)
    }else if(input$input_mode == "Built-in anatomies"){
      req(input$choice)
      if (input$choice == "root cross section"){
        section = get_root_section("./www/Arabido4bis.xml")%>%
          mutate(x = x*0.16, y= y*0.16)
      }
      if (input$choice == "root tip"){
        section = get_root_section("./www/root_tip.xml")%>%
          mutate(x = x*0.16, y= y*0.16)
      }
      if(input$choice == "apical meristem"){
        section = get_root_section("./www/shoot_apex.xml")%>%
          mutate(x = x/2, y= y/2)
      }
      if(input$choice == "embryon"){
        section = get_root_section("./www/embryon.xml")%>%
          mutate(x = x*0.16, y= y*0.16)
      }
      if(input$choice == "pavement cells"){
        section = NULL
        roi_sampl = list.files("./www/",pattern = ".roi", full.names = T)
        for (fl in roi_sampl) {
            section <- rbind(section, get_cell(fl))
        }
        binder = tibble(id_cell = sort(unique(section$id_cell)), 
                        type = c("pavement", "pavement", "pavement",
                                 "C1", "M", "C1", "C2", "M", "pavement", "pavement", "pavement",
                                 "M", "GC", "GC"))
        section = section%>%mutate(x =  x/2.5, y = -y/2.5)%>%
          left_join(binder, by = "id_cell")
        id_global = tibble(id_global  = 1:length(unique(section$id_cell)), id_cell = unique(section$id_cell))
        section = left_join(section, id_global, by = "id_cell")%>%
          mutate(id_cell = id_global)%>%select(-id_global)
      }
      
    }else if(input$input_mode == "Upload ROIs"){
      req(input$all_roi_files$datapath)
      section = NULL
      for(fl in input$all_roi_files$datapath){
        section = rbind(section, get_cell(fl))
      }
      section = section%>%mutate(type = "general", x =  -x/input$dpi, y = -y/input$dpi)
    }
    return(section)
  })
  
  # Show original anatomy
  output$anatomy_plot <- renderPlot({
    req(root_data())
    print("anatomical data loaded")
    root_data()%>%
      ggplot()+
      geom_polygon(aes(x,y, group = id_cell, fill = type), colour = "white")+
      viridis::scale_fill_viridis(discrete = T)+
      theme_classic()+
      coord_fixed()+
      theme(axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank())
  })
  
  # Run geometry prep
  geo_data <- eventReactive(input$prep_btn, {
    req(root_data(),input$input_mode, input$choice,input$thickness, input$corner)
    section = root_data()
    if(input$input_mode == "Upload XML"){
      geo = prep_geo(section, cell_wall_thickness = input$thickness, corner_smoothing = input$corner)
    }else if(input$input_mode == "Upload ROIs"){
      geo = prep_geo(section, cell_wall_thickness = 
                   tibble(type =c("outerwall", "default"),
                          value = c(input$thickness+1,input$thickness)), 
                 corner_smoothing = input$corner)
      }else{
        geo = prep_geo(section, cell_wall_thickness = 
                   tibble(type =c("outerwall", "default"),
                          value = c(input$thickness,input$thickness)), 
                 corner_smoothing = input$corner)
      }
    return(geo)
  })
  

  # Show preprocessed geometry
  output$geo_plot <- renderPlot({
    req(geo_data())
    
    agr = geo_data()%>%filter(id_cell == max(id_cell))
    
    scale_bar = c(floor(quantile(agr$x)[4]/10)*10, floor(max(agr$x)/10)*10, min(agr$y)-sd(agr$y)/10, min(agr$y)-2*sd(agr$y)/10)
    scale_dist = scale_bar[2]-scale_bar[1]
    
    center = c(mean(geo_data()$x), mean(geo_data()$y))
    central = geo_data() %>% filter(id_cell != max(id_cell)) %>% 
      dplyr::group_by(id_cell) %>% 
      dplyr::summarise(mx = mean(x), my = mean(y), .groups = "drop")%>%
      mutate(euc = sqrt((mx-center[1])^2+(my-center[2])^2))
    center_cell_id = central$id_cell[central$euc == min(central$euc)][1]
    
    ggplot()+
      geom_polygon(aes(x,y, group = id_cell), size = 1,colour = "darkgreen", data = geo_data()%>%filter(id_cell == max(id_cell)))+
      geom_polygon(aes(x,y, group = id_cell), size = 0.1,colour = "red", fill = "white", data = geo_data()%>%filter(id_cell != max(id_cell),
                                                                                                                    id_cell != center_cell_id ))+
      geom_polygon(aes(x,y, group = id_cell), size = 0.2,colour = "blue", fill = "white", data = geo_data()%>%filter(id_cell == center_cell_id ))+
      viridis::scale_fill_viridis()+
      geom_segment(aes(x = scale_bar[1], xend = scale_bar[2], y = scale_bar[3], yend = scale_bar[3]),
                   linewidth = 2)+
      geom_text(aes(x = (scale_bar[2]+scale_bar[1])/2, y= scale_bar[4]),
                    label = paste0(scale_dist))+
      theme_classic()+
      coord_fixed()+
      theme(axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank())
  })
  # 
  # # Run geometry prep
  # geo_data <- eventReactive(input$prep_btn, {
  #   
  # })
  
  # Export logic
  output$download_btn <- downloadHandler(
    filename = function() {
      switch(input$export_type,
             "PNG" = "plot.png",
             "SVG" = "plot.svg",
             "CSV" = "geometry.csv",
             "GEO" = "geometry.geo"
      )
    },
    
    content = function(file) {
      switch(input$export_type,
             "PNG" = {
               png(file, width = 800, height = 800)
               root_data()%>%
                 ggplot()+
                 geom_polygon(aes(x,y, group = id_cell, fill = type), colour = "white", alpha = 0.5)+
                 viridis::scale_fill_viridis(discrete = T)+
                 theme_classic()+
                 coord_fixed()+
                 theme(axis.line=element_blank(),
                       axis.text.x=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks=element_blank(),
                       axis.title.x=element_blank(),
                       axis.title.y=element_blank())
               dev.off()
             },
             
             "SVG" = {
               svg(file, width = 8, height = 8)
               root_data()%>%
                 ggplot()+
                 geom_polygon(aes(x,y, group = id_cell, fill = type), colour = "white", alpha = 0.5)+
                 viridis::scale_fill_viridis(discrete = T)+
                 theme_classic()+
                 coord_fixed()+
                 theme(axis.line=element_blank(),
                       axis.text.x=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks=element_blank(),
                       axis.title.x=element_blank(),
                       axis.title.y=element_blank())
               dev.off()
             },
             
             "CSV" = {
               utils::write.csv(geo_data(), file, row.names = FALSE)
             },
             
             "GEO" = {
               write_geo(data = geo_data(), dim = 2, path_geo = file, Celldomain =F)
             }
      )
    }
  )
}
# Run the application 
shinyApp(ui = ui, server = server)


