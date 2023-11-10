BEGIN {
    OFS="\t"
    padding=50
}
{
    if ($6 == "+") {
        left="D";
        right="A";
    } else {
        left="A";
        right="D";
    }
    print $1,$2-padding,$2+padding,left,$4,$5,$6,$7;
    print $1,$3-padding,$3+padding,right,$4,$5,$6,$7
}
