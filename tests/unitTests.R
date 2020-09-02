# Our package. Used for the test suite name
pkgname <- c("MOGAMUN", "RUnit")

for (p in pkgname) {
    require(p, quietly=TRUE, character.only=TRUE) || 
    stop("package '", p, "' not found")
}

pkgname <- "MOGAMUN"

# Determine which files to load (have to start with test_ and end with .R)
pattern <- "^test.*\\.R$"

# Which functions to run. Have to start with 'test.'
testFunctionRegexp = "^test.+"

# Path to the unit tests folder in the package
directory <- paste0(system.file("tests", package=pkgname), "/")
if (.Platform$OS.type == "windows") { # windows needs the absolute file names
    directory <- list.files(directory, pattern = pattern, full.names = TRUE)
} 

for (dir in directory) {
    # Define RUnit test suite
    suite <- defineTestSuite(name=paste(pkgname, "RUnit Tests"),
                             dirs=dir,
                             testFileRegexp=pattern,
                             testFuncRegexp = testFunctionRegexp,
                             rngKind="default",
                             rngNormalKind="default")
    
    # Run tests
    result <- runTestSuite(suite)
    
    # Display result tests on the console
    printTextProtocol(result)
    
    # Write results in JUnit-like xml format
    printJUnitProtocol(result, fileName="junit.xml")
}

