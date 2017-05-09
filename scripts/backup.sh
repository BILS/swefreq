#!/bin/sh

# This script is running hourly on Hertz to backup the user database and
# the downloadable tar archives.
# The cronjob belongs to user "andkaha" and is simply
#   @hourly /data/SweFreq/backup.sh >/dev/null 2>&1
#
# The script may also be run manually on the command line.

PATH="/bin:/usr/bin:$PATH"
export PATH

data_home="/data/SweFreq"

userdb_base="$data_home/userdb-backup"
userdb_dir="$userdb_base/$( date '+%Y-%m' )"
userdb_file="$userdb_dir/tornado-userdb.$( date '+%Y%m%d-%H%M%S' ).dump"

release_backups="$data_home/data-backup/release"
container_dir="/var/lib/lxd/containers/swefreq-web/rootfs/opt/swefreq/tornado/release"

if [ ! -d "$userdb_base" ]; then
    printf '"%s" missing\n' "$userdb_base" >&2
    exit 1
elif [ ! -d "$release_backups" ]; then
    printf '"%s" missing\n' "$release_backups" >&2
    exit 1
fi

trap 'rm -f "$tmpbackup"' EXIT
tmpbackup="$( mktemp -p "$userdb_base" )"

# Dump database, and remove the "Dump completed" comment at the end to
# be able to compare with previous dump.
lxc exec swefreq-web -- mysqldump -u swefreq swefreq |
sed '/^-- Dump completed on/d' >"$tmpbackup"

gzip --best "$tmpbackup"

if [ ! -d "$userdb_dir" ] ||
   [ ! -e "$userdb_base/latest.dump.gz" ] ||
   ! zcmp -s "$tmpbackup.gz" "$userdb_base/latest.dump.gz"
then
    mkdir -p "$userdb_dir"
    mv "$tmpbackup.gz" "$userdb_file.gz"

    # Create a symbolic link to the latest backup
    ln -sf "$userdb_file.gz" "$userdb_base/latest.dump.gz"

    echo 'New userdb backup made'
else
    echo 'No new userdb backup needed'
fi

# Use rsync to sync the "release" directory
rsync --archive --verbose --progress "$container_dir/" "$release_backups/"
