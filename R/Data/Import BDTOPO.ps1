# Choose output Directory
$outputdir = 'Documents/Stage Vincent/Test Download/'

# url from which to download (BDTOPO)
$url       = 'https://geoservices.ign.fr/bdtopo'

# Selected year for data
$year = '2008'

# Index of the file to start download
$startFile = 7

# Create needed directories
$directoriesNeeded = @(-join('C:\Users\vincent\Documents\Stage Vincent\Test download\',$year,'_Full'),-join('C:\Users\vincent\Documents\Stage Vincent\Test download\',$year,'_Compressed_Files'),-join('C:\Users\vincent\Documents\Stage Vincent\Test download\',$year,'_Batiment_Only'))
foreach($directoryToCreate in $directoriesneeded){
    if (-not (Test-Path -LiteralPath $directoryToCreate)) {
        
        try {
            New-Item -Path $directoryToCreate -ItemType Directory -ErrorAction Stop | Out-Null #-Force
        }
        catch {
            Write-Error -Message "Unable to create directory " $directoryToCreate ". Error was: $_" -ErrorAction Stop
        }
    }
}

# Enable TLS 1.2 and TLS 1.1 protocols
[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12, [Net.SecurityProtocolType]::Tls11

# Get all available urls in the website
$WebResponse = Invoke-WebRequest -Uri $url

# Get the list of links for Departements of the selected year
$Files_to_download = $WebResponse.Links | Where-Object {$_.innerHTML -like -join('*_',$year,'-*') -and $_.innerHTML -like'*D0*'} | Select-Object -ExpandProperty href -Skip ($startFile-1)

# Start downloading and unzipping
$compt = $startFile
foreach($7zfile in $Files_to_download){
    Write-Host "Downloading file "$compt.ToString()"/"$Files_to_download.Count.ToString()
    Write-Host $7zfile
    $filepath = -join($outputdir,$year,"_Compressed_Files/","Departement_",$compt.ToString(),"_",$year,".7z")
    Invoke-WebRequest -Uri $7zfile -OutFile $filepath

    Write-Host "Unzipping file..."
    $destinationUnzipPath = -join($outputdir,$year,"_Full/","Departement_",$compt.ToString(),"_",$year)
    
    # Unzip using 7z.exe
    & ${env:ProgramFiles}\7-Zip\7z.exe x $filepath "-o$($destinationUnzipPath)" -y > $null

    # Get the BATIMENT.shp file and copy it in a specific directory
    $finaldestination = -join($outputdir,$year,"_Batiment_Only/")
    $batimentfiles = Get-ChildItem $destinationUnzipPath -Force -Recurse | Where-Object{$_.Name -like '*BATI*' -and ($_.Name -like '*.shp' -or $_.Name -like '*.cpg' -or $_.Name -like '*.dbf' -or $_.Name -like '*.prj' -or $_.Name -like '*.shx')}
    for (($i=0); $i -lt $batimentfiles.Count; $i++){
        $batfile_fullpath = $batimentfiles.FullName[$i]
        $batshpfile = $batimentfiles.Name[$i]
        Copy-Item $batfile_fullpath -Destination $finaldestination
        Rename-Item -Path (-join($finaldestination,$batshpfile)) -NewName (-join("Departement_",$compt.ToString(),"_",$year,"_",$batshpfile))
    }
    $compt++
}


