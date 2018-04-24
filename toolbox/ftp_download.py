import sys
import os
import ftplib

def main(argv):
    ftp_ip_path = argv[1]
    ftp_user = argv[2]
    ftp_pass = argv[3]
    #secret_key = argv[4]
    sep = ftp_ip_path.split(':')
    ftp_ip = sep[0]
    ftp_path = sep[1]
    #decrypt_cmd = "openssl enc -d -aes-256-cbc -in %s -out %s -k "+secret_key
    ftp = ftplib.FTP(ftp_ip)
    ftp.login(ftp_user,ftp_pass)
    ftp.cwd(ftp_path)
    files = []
    ftp.dir(files.append)
    files = [f.split(' ')[-1] for f in files]
    downloaded = os.listdir('.')
    down_count = 0
    up_count = 0
    for f in files:
        if f in downloaded:
            print "already downloaded",f
            down_count += 1
            continue
        print "downloading",f
        up_count += 1
        with open(f,'wb') as local_file:
            ftp.retrbinary('RETR %s' % f, local_file.write)
        #cmd = decrypt_cmd % (f,f.split('.enc')[0])
        #cmd = cmd.split(' ')
        #try:
        #    subprocess.check_call(cmd)
        #except:
        #    print "failed to decrypt",f 
    try:
        ftp.quit()
    except:
        ftp.close()
    print "Skipped:",down_count
    print "Downloaded:",up_count


if __name__ == '__main__':
    main(sys.argv)

