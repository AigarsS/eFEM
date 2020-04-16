
' Sub sijaGEM()

'     skaitsL = Range("A8").Value
'     skaitsS = Cells(11 + skaitsL, 1).Value
'     I = Cells(19 + skaitsL + skaitsS, 13).Value
'     E = Cells(20 + skaitsL + skaitsS, 13).Value
    
'     ReDim l_e(1 To skaitsL * 10 + 1)
'     ReDim F_x(1 To skaitsS)
    
'     'Nolasa x koordinātu punktiem, kas sadala pa laidumiem
'     l_e(1) = 0
'     For m = 1 To skaitsL
'         l_e(10 * m - 8) = 1 * Cells(8 + m, 13).Value / 10 + l_e(10 * m - 9)
'         l_e(10 * m - 7) = 2 * Cells(8 + m, 13).Value / 10 + l_e(10 * m - 9)
'         l_e(10 * m - 6) = 3 * Cells(8 + m, 13).Value / 10 + l_e(10 * m - 9)
'         l_e(10 * m - 5) = 4 * Cells(8 + m, 13).Value / 10 + l_e(10 * m - 9)
'         l_e(10 * m - 4) = 5 * Cells(8 + m, 13).Value / 10 + l_e(10 * m - 9)
'         l_e(10 * m - 3) = 6 * Cells(8 + m, 13).Value / 10 + l_e(10 * m - 9)
'         l_e(10 * m - 2) = 7 * Cells(8 + m, 13).Value / 10 + l_e(10 * m - 9)
'         l_e(10 * m - 1) = 8 * Cells(8 + m, 13).Value / 10 + l_e(10 * m - 9)
'         l_e(10 * m) = 9 * Cells(8 + m, 13).Value / 10 + l_e(10 * m - 9)
'         l_e(10 * m + 1) = Cells(8 + m, 13).Value + l_e(10 * m - 9)
'     Next m

'     'Nolasa x koordinātu punktiem, kuros pieliktas slodzes
'     n = 0
'     For m = 1 To skaitsS
'         slodzesTips = Cells(12 + skaitsL + m, 4).Value
'         Select Case slodzesTips
'             Case "Punktveida"
'                 F_x(m) = Cells(12 + skaitsL + m, 10).Value
'             Case "Linijveida"
'                 F_x(m) = Cells(12 + skaitsL + m, 10).Value
'                 n = n + 1
'                 ReDim Preserve F_x(1 To skaitsS + n)
'                 F_x(skaitsS + n) = Cells(12 + skaitsL + m, 13).Value
'         End Select
'     Next m

'     'saliek visus punktus vienā matricā
'     ReDim xtemp(1 To UBound(F_x, 1) + UBound(l_e, 1))
'     ReDim x(1 To UBound(F_x, 1) + UBound(l_e, 1))
'     n = 0
'     For m = 1 To UBound(l_e, 1)
'         xtemp(m) = l_e(m)
'     Next m
'     For m = UBound(l_e, 1) + 1 To UBound(xtemp, 1)
'         n = n + 1
'         xtemp(m) = F_x(n)
'     Next m
    
'     'sašķiro matricas elementus no mazākā uz lielāko
'     Call QuickSort(xtemp, LBound(xtemp), UBound(xtemp))
    
'     'izņem no matricas vienādos elementus
'     n = 0
'     For m = 1 To UBound(xtemp, 1) - 1
'         If m = UBound(xtemp, 1) - 1 And xtemp(m) <> xtemp(m + 1) Then
'             n = n + 1
'             x(n) = xtemp(m)
'             x(n + 1) = xtemp(m + 1)
'         ElseIf m = UBound(xtemp, 1) - 1 And xtemp(m) = xtemp(m + 1) Then
'             x(n + 1) = xtemp(m)
'         ElseIf xtemp(m) <> xtemp(m + 1) Then
'             n = n + 1
'             x(n) = xtemp(m)
'         End If
'     Next m
'     ReDim Preserve x(1 To n + 1)

'     'izveido globālo stinguma matricu
'     a = 0
'     ReDim K(1 To 2 * UBound(x, 1), 1 To 2 * UBound(x, 1))
'     ReDim stiffK(1 To UBound(K, 1), 1 To UBound(K, 2))
'     For m = 1 To UBound(x, 1) - 1
'         L = x(m + 1) - x(m)
'         k_e = EulerBernoulli(L)
'         S = E * I / L ^ 3
'         For n = 1 To 4
'             For j = 1 To 4
'               K(n + a, j + a) = K(n + a, j + a) + S * k_e(n, j)
'               stiffK(n + a, j + a) = stiffK(n + a, j + a) + S * k_e(n, j)
'             Next j
'         Next n
'         a = a + 2
'     Next m
    
'     'Robežnosacījumu pielikšana
'     ReDim FBcs(1 To UBound(K, 1), 1 To 1)
'     For a = 0 To skaitsL
'         balstijums = Cells(9 + a, 10).Value
'         If a = 0 Then
'             x1 = 0
'         Else: x1 = x1 + Cells(8 + a, 13).Value
'         End If
'         Select Case balstijums
'             Case "Brivi balstits"
'                 For m = 1 To UBound(x, 1)
'                     If x1 = x(m) Then
'                         For n = 1 To UBound(K, 1)
'                             K(2 * m - 1, n) = 0
'                         Next n
'                         K(2 * m - 1, 2 * m - 1) = 1
'                         FBcs(2 * m - 1, 1) = "BCS"
'                     End If
'                 Next m
'             Case "Iespilejums"
'                 For m = 1 To UBound(x, 1)
'                     If x1 = x(m) Then
'                         For n = 1 To UBound(K, 1)
'                             K(2 * m - 1, n) = 0
'                             K(2 * m, n) = 0
'                         Next n
'                          K(2 * m - 1, 2 * m - 1) = 1
'                          K(2 * m, 2 * m) = 1
'                          FBcs(2 * m - 1, 1) = "BCS"
'                          FBcs(2 * m, 1) = "BCS"
'                     End If
'                 Next m
'         End Select
'     Next a

'     'Tukšajās vietās ievieto "0"
'     For m = 1 To UBound(K, 1)
'         For n = 1 To UBound(K, 2)
'             If IsEmpty(K(m, n)) Then K(m, n) = 0
'             If IsEmpty(stiffK(m, n)) Then stiffK(m, n) = 0
'         Next n
'     Next m
      
    'Inversā matrica un atrisinājums ar pārvietojumu matricu
    invK = WorksheetFunction.MInverse(K)

    'Slodžu matrica
    ReDim d(1 To skaitsS, 1 To UBound(K, 1))
    ReDim F(1 To UBound(K, 1), 1 To 1)
    ReDim F0(1 To skaitsS, 1 To UBound(K, 1))
    ReDim M0(1 To skaitsS, 1 To UBound(K, 1))
    For a = 1 To skaitsS
        slodzesTips = Cells(12 + skaitsL + a, 4).Value
        For m = 1 To UBound(F, 1)
            F(m, 1) = 0
        Next m
        Select Case slodzesTips
            Case "Punktveida"
                tempF = Cells(12 + skaitsL + a, 7).Value
                x1 = Cells(12 + skaitsL + a, 10).Value
                For m = 1 To UBound(x, 1)
                    If x1 = x(m) Then F(m * 2 - 1, 1) = F(m * 2 - 1, 1) - tempF
                Next m
                
                For m = 1 To UBound(x, 1) - 1
                    If x(m) >= x1 And x(m) < x2 Then
                        F0(a, m * 2 - 1) = 0
                        F0(a, m * 2) = 0
                        M0(a, m * 2 - 1) = 0
                        M0(a, m * 2) = 0
                    End If
                Next m
            Case "Linijveida"
                tempF = Cells(12 + skaitsL + a, 7).Value
                x1 = Cells(12 + skaitsL + a, 10).Value
                x2 = Cells(12 + skaitsL + a, 13).Value
                For m = 1 To UBound(x, 1) - 1
                    If x(m) >= x1 And x(m) < x2 Then
                        F(m * 2 - 1, 1) = F(m * 2 - 1, 1) - tempF * (x(m + 1) - x(m)) / 2
                        F(m * 2, 1) = F(m * 2, 1) - tempF * (x(m + 1) - x(m)) ^ 2 / 12
                        F(m * 2 + 1, 1) = F(m * 2 + 1, 1) - tempF * (x(m + 1) - x(m)) / 2
                        F(m * 2 + 2, 1) = F(m * 2 + 2, 1) + tempF * (x(m + 1) - x(m)) ^ 2 / 12
                    End If
                Next m
                
                For m = 1 To UBound(x, 1) - 1
                    If x(m) >= x1 And x(m) < x2 Then
                        F0(a, m * 2 - 1) = -tempF * (x(m + 1) - x(m)) / 2
                        F0(a, m * 2) = -tempF * (x(m + 1) - x(m)) / 2
                        M0(a, m * 2 - 1) = -tempF * (x(m + 1) - x(m)) ^ 2 / 12
                        M0(a, m * 2) = tempF * (x(m + 1) - x(m)) ^ 2 / 12
                    End If
                Next m
        End Select

        For m = 1 To UBound(F, 1)
            If FBcs(m, 1) = "BCS" Then F(m, 1) = 0
        Next m
        
        tempd = WorksheetFunction.MMult(invK, F)
        
        For m = 1 To UBound(tempd, 1)
            d(a, m) = tempd(m, 1)
        Next m
        
    Next a
    
    Kd = WorksheetFunction.MMult(stiffK, tempd)

    
'    pārvietojumi
    mesh = 100
    ReDim u(1 To skaitsS, 1 To 1)
    ReDim u_x(1 To 1)
    a = 1
    For m = 2 To UBound(x, 1)
        L = x(m) - x(m - 1)
        For x1 = 0 To L Step mesh
            N1 = 1 / L ^ 3 * (2 * x1 ^ 3 - 3 * x1 ^ 2 * L + L ^ 3)
            N2 = 1 / L ^ 3 * (x1 ^ 3 * L - 2 * x1 ^ 2 * L ^ 2 + x1 * L ^ 3)
            N3 = 1 / L ^ 3 * (-2 * x1 ^ 3 + 3 * x1 ^ 2 * L)
            N4 = 1 / L ^ 3 * (x1 ^ 3 * L - x1 ^ 2 * L ^ 2)
            For b = 1 To skaitsS
                u(b, a) = N1 * d(b, 2 * m - 3) + N2 * d(b, 2 * m - 2) + N3 * d(b, 2 * m - 1) + N4 * d(b, 2 * m)
            Next b
            If x1 = 0 And m <> 2 Then
                u_x(a) = x(m - 1) - x(1)
            ElseIf x1 = 0 And m = 2 Then
                u_x(a) = x1
            ElseIf x1 = L And m <> UBound(x, 1) Then
                a = a - 1
            ElseIf x1 = L And m = UBound(x, 1) Then
               u_x(a) = u_x(a - 1) + mesh
               Exit For
            Else
                u_x(a) = u_x(a - 1) + mesh
            End If
            a = a + 1
            ReDim Preserve u_x(1 To a)
            ReDim Preserve u(1 To skaitsS, 1 To a)
            If L - x1 < mesh And m = UBound(x, 1) Then
                u_x(a) = u_x(a - 1) + (L - x1)
                x1 = L
                N1 = 1 / L ^ 3 * (2 * x1 ^ 3 - 3 * x1 ^ 2 * L + L ^ 3)
                N2 = 1 / L ^ 3 * (x1 ^ 3 * L - 2 * x1 ^ 2 * L ^ 2 + x1 * L ^ 3)
                N3 = 1 / L ^ 3 * (-2 * x1 ^ 3 + 3 * x1 ^ 2 * L)
                N4 = 1 / L ^ 3 * (x1 ^ 3 * L - x1 ^ 2 * L ^ 2)
                For b = 1 To skaitsS
                    u(b, a) = N1 * d(b, 2 * m - 3) + N2 * d(b, 2 * m - 2) + N3 * d(b, 2 * m - 1) + N4 * d(b, 2 * m)
                Next b
            End If
        Next x1
    Next m
    
    'Šķērsspēki
    ReDim V_x(1 To skaitsS, 1 To (UBound(x, 1) - 1) * 2)
    ReDim M_x(1 To skaitsS, 1 To (UBound(x, 1) - 1) * 2)
    ReDim koordVx(1 To UBound(x, 1) * 2, 1 To 1)
    koordVx(1, 1) = 0
    n = 1
    For m = 1 To UBound(x, 1) - 1
        L = x(m + 1) - x(m)
        k_e = EulerBernoulli(L)
        S = E * I / L ^ 3
        If m <> 1 Then
            koordVx(n, 1) = koordVx(n - 1, 1)
        End If
        koordVx(n + 1, 1) = koordVx(n, 1) + L
        
        For b = 1 To skaitsS

            V_x(b, n) = V_x(b, n) + S * (k_e(1, 1) * d(b, 2 * m - 1) + k_e(1, 2) * d(b, 2 * m) _
            + k_e(1, 3) * d(b, 2 * m + 1) + k_e(1, 4) * d(b, 2 * m + 2)) - F0(b, 2 * m - 1)
    
            V_x(b, n + 1) = V_x(b, n + 1) - S * (k_e(3, 1) * d(b, 2 * m - 1) + k_e(3, 2) * d(b, 2 * m) _
            + k_e(3, 3) * d(b, 2 * m + 1) + k_e(3, 4) * d(b, 2 * m + 2)) + F0(b, 2 * m)
    
            M_x(b, n) = M_x(b, n) + S * (k_e(2, 1) * d(b, 2 * m - 1) + k_e(2, 2) * d(b, 2 * m) _
            + k_e(2, 3) * d(b, 2 * m + 1) + k_e(2, 4) * d(b, 2 * m + 2)) - M0(b, 2 * m - 1)
    
            M_x(b, n + 1) = M_x(b, n + 1) - S * (k_e(4, 1) * d(b, 2 * m - 1) + k_e(4, 2) * d(b, 2 * m) _
            + k_e(4, 3) * d(b, 2 * m + 1) + k_e(4, 4) * d(b, 2 * m + 2)) + M0(b, 2 * m)
        
        Next b

        n = n + 2
    Next m
'
    ' Worksheets(2).UsedRange.ClearContents
'    For m = 1 To UBound(u, 1)
'        Worksheets(2).Cells(5 + m, 3).Value = u_x(m)
'        Worksheets(2).Cells(5 + m, 5).Value = u(m)
'    Next m
'
'    For m = 1 To UBound(Kd, 1)
'        Worksheets(2).Cells(5 + m, 8).Value = Kd(m, 1)
'    Next m
'
'    For m = 1 To UBound(d, 1)
'        For n = 1 To UBound(d, 2)
'            Worksheets(2).Cells(5 + n, 10 + m).Value = d(m, n)
'        Next n
'    Next m
    
'     For m = 1 To UBound(u_x, 1)
'         Worksheets(2).Cells(m, 1).Value = u_x(m)
'     Next m
    
'     For m = 1 To skaitsS
'         For n = 1 To UBound(u, 2)
'             Worksheets(2).Cells(n, 2).Value = Worksheets(2).Cells(n, 2).Value + u(m, n)
'         Next n
'     Next m
    
'     For m = 1 To UBound(koordVx, 1)
'         Worksheets(2).Cells(m, 5).Value = koordVx(m, 1)
'     Next m

'     For m = 1 To skaitsS
'         For n = 1 To UBound(V_x, 2)
'             Worksheets(2).Cells(n, 6).Value = Worksheets(2).Cells(n, 6).Value + V_x(m, n) / 10 ^ 3
'             Worksheets(2).Cells(n, 7).Value = Worksheets(2).Cells(n, 7).Value + M_x(m, n) / 10 ^ 6
'         Next n
'     Next m
' '
'    For m = 1 To UBound(V_x, 1)
'        Worksheets(2).Cells(5 + m, 12).Value = koordVx(m, 1)
'        Worksheets(2).Cells(5 + m, 13).Value = V_x(m) * 10 ^ -3
'        Worksheets(2).Cells(5 + m, 14).Value = M_x(m) * 10 ^ -6
'    Next m
'
'    For m = 1 To UBound(x, 1)
'        Worksheets(2).Cells(5 + m, 15).Value = x(m)
'    Next m

End Sub

' Sub QuickSort(arr, Lo As Long, Hi As Long)
'   Dim varPivot As Variant
'   Dim varTmp As Variant
'   Dim tmpLow As Long
'   Dim tmpHi As Long
'   tmpLow = Lo
'   tmpHi = Hi
'   varPivot = arr((Lo + Hi) \ 2)
'   Do While tmpLow <= tmpHi
'     Do While arr(tmpLow) < varPivot And tmpLow < Hi
'       tmpLow = tmpLow + 1
'     Loop
'     Do While varPivot < arr(tmpHi) And tmpHi > Lo
'       tmpHi = tmpHi - 1
'     Loop
'     If tmpLow <= tmpHi Then
'       varTmp = arr(tmpLow)
'       arr(tmpLow) = arr(tmpHi)
'       arr(tmpHi) = varTmp
'       tmpLow = tmpLow + 1
'       tmpHi = tmpHi - 1
'     End If
'   Loop
'   If Lo < tmpHi Then QuickSort arr, Lo, tmpHi
'   If tmpLow < Hi Then QuickSort arr, tmpLow, Hi
' End Sub

' Function EulerBernoulli(Le) As Variant
'     Dim kLoc(1 To 4, 1 To 4) As Double
'     kLoc(1, 1) = 12:         kLoc(1, 2) = 6 * Le:        kLoc(1, 3) = -12:         kLoc(1, 4) = 6 * Le
'     kLoc(2, 1) = 6 * Le:     kLoc(2, 2) = 4 * Le ^ 2:    kLoc(2, 3) = -6 * Le:     kLoc(2, 4) = 2 * Le ^ 2
'     kLoc(3, 1) = -12:        kLoc(3, 2) = -6 * Le:       kLoc(3, 3) = 12:          kLoc(3, 4) = -6 * Le
'     kLoc(4, 1) = 6 * Le:     kLoc(4, 2) = 2 * Le ^ 2:    kLoc(4, 3) = -6 * Le:     kLoc(4, 4) = 4 * Le ^ 2
'     EulerBernoulli = kLoc
' End Function