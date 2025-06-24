pipeline {
    agent any

    environment {
        R_LIBS_USER = "${WORKSPACE}/Rlibs"  // custom R lib path
    }

    stages {
        stage('Clone Repo') {
            steps {
                git 'https://github.com/jcaperella29/RNA_SEQ_APP.git'
            }
        }

        stage('Install R Packages') {
            steps {
                sh '''
                    mkdir -p $R_LIBS_USER
                    Rscript -e "install.packages(c('shiny', 'shinydashboard', 'DT', 'ggplot2', 'plotly', 'dplyr', 'lintr'), repos='https://cran.rstudio.com')"
                '''
            }
        }

        stage('Lint R Code') {
            steps {
                sh '''
                    Rscript -e "lintr::lint_dir()"
                '''
            }
        }

        // Optional: Add testing or deployment stages here
    }
}
