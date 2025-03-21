

library(blastula)

log_file <- "error_log.txt"

# Only send if log exists and is not empty
if (file.exists(log_file) && file.info(log_file)$size > 0) {
  
  # Read the log contents
  log_text <- paste(readLines(log_file), collapse = "\n")
  
  # Compose the email
  email <- compose_email(
    body = md(paste0("
## 🛠️ Daily Shiny App Error Log

The following errors were logged today:
"))
  )
  
  # Send the email
  smtp_send(
    email,
    from = "jcaperella@gmail.com",
    to = "jcaperella@gmail.com",
    subject = paste("Shiny Error Log -", Sys.Date()),
    credentials = creds_key(id = "gmail_creds")
  )
  
  # Archive the log
  dir.create("logs", showWarnings = FALSE)
  file.rename(log_file, paste0("logs/log_", Sys.Date(), ".txt"))
}
