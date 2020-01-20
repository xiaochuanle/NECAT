count=$1
action=$2
cfgfile=$3

necat_path=$(cd "$(dirname "$0")";pwd)

for ((i=0;i<$count;i++)) do
    $necat_path/necat.pl $action $cfgfile
    temp=$?
    echo "return $temp"
    if [ $temp -eq 0 ]; then
        break
    fi
done

