from curses import meta
import os 
import argparse
import subprocess
import pandas as pd 

parser = argparse.ArgumentParser(description='NCBI download reference genome for taxid of tax name',
                                 prog='NCBI_download_ref',
                                 usage='')
parser.version = '1.2'
parser.add_argument('--name', help= "Organism Name",type=str)  
parser.add_argument('--txid', help= "Taxonomy ID",type=str)  
parser.add_argument('--id', help= "Taxonomy ID",type=str)  
parser.add_argument('--kraken-report', help= "Kraken report path",type=str)  
parser.add_argument("-o","--outdir",help= "Output name ( output dir )")
parser.add_argument("--direct",help= "No user choise, direct download first record", action='store_true')

args = parser.parse_args()

org_name = ""
output = ""
directdownload = args.direct
def retrive_ref(tax,taxid=False):
    print(f"\033[;32;1mSearching for {tax}...\033[;39;m")
    if taxid == False:
        out = subprocess.Popen(f"esearch -db assembly -query '{tax}[ORGN]' | esummary | xtract -pattern DocumentSummary \
                -def NA -element FtpPath_GenBank,Organism,AssemblyAccession,Taxid,assembly-status -block Stat \
                -if Stat@category -equals contig_count -or Stat@category -equals contig_l50 \
                -or Stat@category -equals contig_n50 -or Stat@category -equals total_length \
                -element Stat", shell=True, stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
    else:
        out = subprocess.Popen(f"esearch -db genome -query '{tax}' |elink -target assembly|esummary| xtract -pattern DocumentSummary \
                -def NA -element FtpPath_GenBank,Organism,AssemblyAccession,Taxid,assembly-status -block Stat \
                -if Stat@category -equals contig_count -or Stat@category -equals contig_l50 \
                -or Stat@category -equals contig_n50 -or Stat@category -equals total_length \
                -element Stat", shell=True, stdout=subprocess.PIPE).communicate()[0].decode('utf-8')

    raw_lst = []
    lst = out.strip().split("\n")
    for i in lst:
        raw_lst.append(i.split("\t"))
    df = pd.DataFrame(raw_lst)
    df.columns =["FtpPath_GenBank",'Organism', 'AssemblyAccession', 
                    'Taxid', 'assembly-status', 'contig_count',
                    'contig_l50', 'contig_n50', 'total_length']

    df["Organism"].unique()
    df = df.astype({"contig_count": int})
    return df


def filter_taxa_df(df,user_input=True):
    global org_name, output, directdownload
    uniques = df["Organism"].unique()
    if user_input is True and not directdownload:
        print("\033[;33;1mPlease Select an Organism: \033[;39;m")
        counter = 1
        for i in uniques:
            print(f"{counter}: {i}")
            counter += 1
        n = int(input("\n"))

    else:
        n = 1
    m = n-1
    out = df[df["Organism"] == df["Organism"].unique()[m]].sort_values("contig_count")
    org_name = str(uniques[m]).replace(" ", "_")
    org_name = org_name.replace("(", "-")
    org_name = org_name.replace(")", "-")
    output = str(args.outdir) + f"/{org_name}.fna.gz"
    uniques = out["assembly-status"].unique()
    if user_input is True and not directdownload:
        print("\033[;33;1mPlease Select an Assembly level: \033[;39;m")
        counter = 1
        for i in uniques:
            print(f"{counter}: {i}")
            counter += 1
        n = int(input("\n"))

    else:
        n = 1
    m = n-1
    out2 = out[out["assembly-status"] == out["assembly-status"].unique()[m]].sort_values("contig_count")
    return out2

def download_ref(df,FILENAME,user_input=True):
    global directdownload
    df2 = df[["FtpPath_GenBank", "AssemblyAccession", "total_length", "assembly-status", "contig_count", "contig_l50", "contig_n50"]]
    df2.index = df2.reset_index(drop = True).index + 1

    if len(df2.index) > 1:
        if user_input is True and not directdownload:
            print("\033[;33;1mPlease Select a Number to Download: \033[;39;m")
            print(df2.loc[:, df2.columns != 'FtpPath_GenBank'].to_string())
            n = int(input("\n"))

        else:
            n = 1
            print("\n" + "\033[;33;1mThe following will be downloaded:\n\033[;37;1m" + str(df2.loc[:1, df2.columns != 'FtpPath_GenBank'].to_string())+ "\033[;39;m\n\n")
        m = n-1
        acc = df2.iloc[m]["AssemblyAccession"]
        ftp = df2.iloc[m]["FtpPath_GenBank"]
    else:
        print(df2.loc[:, df2.columns != 'FtpPath_GenBank'].to_string())
        acc = df2.iloc[0]["AssemblyAccession"]
        ftp = df2.iloc[0]["FtpPath_GenBank"]
    ftp_full = ftp + "/" + ftp.split("/")[-1] + "_genomic.fna.gz"
    os.system(f"wget -cO - {ftp_full} > {FILENAME}")
    print(f"\033[;32;1mDownloaded {acc} @ {FILENAME}\033[;39;m")

def process_kraken(path):
    df = pd.read_table(path)
    df.columns =['perc', 'clade', 'assigned', 'rank', 'txid', 'name']
    out = df[df["rank"] == "S" ].sort_values("perc",ascending=False)
    out.index = out.reset_index(drop = True).index + 1
    name = str(out.iloc[0]["name"]).strip()
    return name

try:
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
except OSError as error:
    print(f"\033[;31;1mCreating output at{args.outdir}: \nraise an error: \033[;39;m{error}")
    exit(1)

if args.name is None :
    if args.txid is None:
        if args.id is None:
            if args.kraken_report is None:
                print("\033[;31;1mNo --name or --txid or --id or --kraken-report supplied: \033[;39;m")
                parser.print_help()
                exit(1)
            else:
                org_name = str(process_kraken(args.kraken_report))
                org = str(process_kraken(args.kraken_report)).replace(" ","_")
                print(f"\033[;33;1mTop species in the supplied report: {org_name}\033[;39;m")
                df = retrive_ref(org_name)
                df_filterred = filter_taxa_df(df,user_input=False)
                download_ref(df_filterred,str(f"{output}/{org}.fna.gz"),user_input=False)
        else:
            df = retrive_ref(args.id,taxid=True)
            df_filterred = filter_taxa_df(df)
            download_ref(df_filterred, str(f"{output}"))
    else:
        df = retrive_ref(str(f"txid{args.txid}"),taxid=True)
        df_filterred = filter_taxa_df(df)
        download_ref(df_filterred, str(f"{output}"))
else:
    df = retrive_ref(args.name)
    df_filterred = filter_taxa_df(df)
    download_ref(df_filterred, str(f"{output}"))