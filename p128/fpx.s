.intel_syntax noprefix

.extern pbits
.set pbytes,200
.set plimbs,25

.section .text

.reduce_once:
    push rbp
    mov rbp, rdi

    mov rdi, [rbp + 0]
    sub rdi, [rip + p + 0]
    mov rsi, [rbp + 8]
    sbb rsi, [rip + p + 8]
    mov rdx, [rbp + 16]
    sbb rdx, [rip + p + 16]
    mov rcx, [rbp + 24]
    sbb rcx, [rip + p + 24]
    mov r8, [rbp + 32]
    sbb r8, [rip + p + 32]
    mov r9, [rbp + 40]
    sbb r9, [rip + p + 40]
    mov r10, [rbp + 48]
    sbb r10, [rip + p + 48]
    mov r11, [rbp + 56]
    sbb r11, [rip + p + 56]

    mov rdi, [rbp + 64]
    sbb rdi, [rip + p + 64]
    mov rsi, [rbp + 72]
    sbb rsi, [rip + p + 72]
    mov rdx, [rbp + 80]
    sbb rdx, [rip + p + 80]
    mov rcx, [rbp + 88]
    sbb rcx, [rip + p + 88]
    mov r8, [rbp + 96]
    sbb r8, [rip + p + 96]
    mov r9, [rbp + 104]
    sbb r9, [rip + p + 104]
    mov r10, [rbp + 112]
    sbb r10, [rip + p + 112]
    mov r11, [rbp + 120]
    sbb r11, [rip + p + 120]

    mov rdi, [rbp + 128]
    sbb rdi, [rip + p + 128]
    mov rsi, [rbp + 136]
    sbb rsi, [rip + p + 136]
    mov rdx, [rbp + 144]
    sbb rdx, [rip + p + 144]
    mov rcx, [rbp + 152]
    sbb rcx, [rip + p + 152]
    mov r8, [rbp + 160]
    sbb r8, [rip + p + 160]
    mov r9, [rbp + 168]
    sbb r9, [rip + p + 168]
    mov r10, [rbp + 176]
    sbb r10, [rip + p + 176]
    mov r11, [rbp + 184]
    sbb r11, [rip + p + 184]

    mov rdi, [rbp + 192]
    sbb rdi, [rip + p + 192]
    setnc al
    movzx rax, al
    neg rax

.global fp_add3
fp_add3:
  push rdi
  call uint_custom_add3
  pop rdi

  jmp .reduce_once

.global fp_add2
fp_add2:
    mov rdx, rsi
    jmp fp_add3

.global fp_sub3
fp_sub3:
  push rdi
  call uint_custom_sub3
  pop rdi

  neg rax

  sub rsp, pbytes
  mov rcx, [rip + p + 0]zzzzzzzzzzz
  and rcx, rax
  mov [rsp + 0], rcx
  .set k, 1
  .rept plimbs-1
      mov rcx, [rip + p + 8*k]
      and rcx, rax
      mov [rsp + 8*k], rcx
      .set k, k+1
  .endr

  mov rcx, [rsp + 0]
  add rcx, [rdi + 0]
  mov [rdi + 0], rcx
  .set k, 1
  .rept plimbs-1
      mov rcx, [rsp + 8*k]
      adc rcx, [rdi + 8*k]
      mov [rdi + 8*k], rcx
      .set k, k+1
  .endr

  add rsp, pbytes
  ret

.global fp_sub2
fp_sub2:
    mov rdx, rsi
    jmp fp_sub3

.global fp_sqr1
fp_sqr1:
    mov rdx, rdi
    mov rsi, rdi
    jmp fp_mul3

.global fp_sqr2
fp_sqr2:
    mov rdx, rsi
    mov rsi, rsi
    jmp fp_mul3

.global fp_pow
fp_pow:
    push rbx
    push r12
    push r13

    mov rbx, rsi
    mov r12, rdi

    lea rsi, [rip + fp_1]
    mov rdi, r12
    call fp_set

.macro POWSTEP, k
        mov r13, [rbx + 8*\k]
        mov rcx, 0

0:
        test r13, 1
        jz 1f

        mov rdi, r12
        mov rsi, r12
        mov rdx, rbx
        call fp_mul3

1:
        mov rdi, r12
        mov rsi, r12
        mov rdx, r12
        call fp_mul3

        shr r13
        inc r12
        test rcx, 64
        jz 0b
.endm

    POWSTEP 0
    POWSTEP 1
    POWSTEP 2
    POWSTEP 3
    POWSTEP 4
    POWSTEP 5
    POWSTEP 6
    POWSTEP 7
    POWSTEP 8
    POWSTEP 9
    POWSTEP 10
    POWSTEP 11
    POWSTEP 12
    POWSTEP 13
    POWSTEP 14
    POWSTEP 15
    POWSTEP 16
    POWSTEP 17
    POWSTEP 18
    POWSTEP 19
    POWSTEP 20
    POWSTEP 21
    POWSTEP 22
    POWSTEP 23
    POWSTEP 24

    pop r13
    pop r12
    pop rbx
    ret

.global fp_inv
fp_inv:
    lea rsi, [rip + p_minus_2]
    jmp fp_pow

.global fp_issquare
fp_issquare:
.set k, 24
.rept plimbs
    pushq [rip + p + 8*k]
    .set k, k-1
.endr

    mov rsi, rsp
    mov cl, 1  /* result */

.l0p:

    .set k, 0
    .set v, 1
    .rept plimbs
        cmpq [rdi + 8*k], v
        jnz 1f
        .set k, k+1
        .set v, 0
    .endr
        jmp .l0pbrk
        1:

        testq [rdi], 1
        jz .shift

    .set k, 24
    .rept 25
        mov rax, [rdi + 8*k]
        cmp rax, [rsi + 8*k]
        ja .reduce
        .if k != 0
        jb .recip
        .endif
        .set k, k-1
    .endr

    .recip:
        mov al, [rdi]
        and al, [rsi]
        and al, 3
        cmp al, 3
        sete al
        xor cl, al
        xchg rdi, rsi

    .reduce:
        mov rax, [rsi +  0]
        sub [rdi +  0], rax
    .set k, 1
    .rept 24
        mov rax, [rsi + 8*k]
        sbb [rdi + 8*k], rax
        .set k, k+1
    .endr

    .shift:
        mov eax, [rsi]
        and eax, 24
        popcnt eax, eax
        test eax, 1
        setz al
        xor cl, al
        shrq [rdi + 192]
        rcrq [rdi + 184]
        rcrq [rdi + 176]
        rcrq [rdi + 168]
        rcrq [rdi + 160]
        rcrq [rdi + 152]
        rcrq [rdi + 144]
        rcrq [rdi + 136]
        rcrq [rdi + 128]
        rcrq [rdi + 120]
        rcrq [rdi + 112]
        rcrq [rdi + 104]
        rcrq [rdi + 96]
        rcrq [rdi + 88]
        rcrq [rdi + 80]
        rcrq [rdi + 72]
        rcrq [rdi + 64]
        rcrq [rdi + 56]
        rcrq [rdi + 48]
        rcrq [rdi + 40]
        rcrq [rdi + 32]
        rcrq [rdi + 24]
        rcrq [rdi + 16]
        rcrq [rdi +  8]
        rcrq [rdi +  0]
        jmp .l0p

.l0pbrk:
    movzx eax, cl
    add rsp, 64
    ret

.global fp_random
fp_random:
    lea rsi, [rip + p]
    jmp uint_custom_random

.global fp_eq
fp_eq:
    xor eax, eax
.set k, 0
.rept plimbs
    mov rdx, [rdi + 8*k]
    xor rdx, [rsi + 8*k]
    or rax, rdx
    .set k, k+1
.endr
    test rax, rax
    setz al
    movzx eax, al
    ret

.global fp_set
fp_set:
    push rdi
    call uint_custom_set
    pop rdi
    mov rsi, rdi
    jmp fp_enc

.global fp_mul3
fp_mul3:
  push rbp
  push rbx

  sub rsp, 216
  mov [rsp+208], rdi
  mov rdi, rsi
  mov rsi, rdx

  xor rax, rax
  mov [rsp+0], rax
  mov [rsp+8], rax
  mov [rsp+16], rax
  mov [rsp+24], rax
  mov [rsp+32], rax
  mov [rsp+40], rax
  mov [rsp+48], rax
  mov [rsp+56], rax
  mov [rsp+64], rax
  mov [rsp+72], rax
  mov [rsp+80], rax
  mov [rsp+88], rax
  mov [rsp+96], rax
  mov [rsp+104], rax
  mov [rsp+112], rax
  mov [rsp+120], rax
  mov [rsp+128], rax
  mov [rsp+136], rax
  mov [rsp+144], rax
  mov [rsp+152], rax
  mov [rsp+160], rax
  mov [rsp+168], rax
  mov [rsp+176], rax
  mov [rsp+184], rax
  mov [rsp+192], rax
  mov [rsp+200], rax

.macro MULSTEP, k,I0,I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,I13,I14,I15,I16,I17,I18,I19,I20,I21,I22,I23,I24

    mov r11,[rsp+\I0]
    mov rdx, [rsi +  0]
    mulx rcx, rdx, [rdi + 8*k]
    add rdx, r11
    mulx rcx, rdx, [rip + inv_min_p_mod_r]

    xor rax, rax

    mulx rbx, rax, [rip + p +  0]
    adox r11, rax
    mov [rsp+\I0], r11

    mov r11,[rsp+\I1]
    mulx rcx, rax, [rip + p + 8]
    adcx r11, rcx
    adox r11, rax
    mov [rsp+\I1], r11

    mov r11,[rsp+\I2]
    mulx rbx, rax, [rip + p + 16]
    adcx r11, rbx
    adox r11, rax
    mov [rsp+\I2], r11

    mov r11,[rsp+\I3]
    mulx rcx, rax, [rip + p + 24]
    adcx r11, rcx
    adox r11, rax
    mov [rsp+\I3], r11

    mov r11,[rsp+\I4]
    mulx rbx, rax, [rip + p + 32]
    adcx r11, rbx
    adox r11, rax
    mov [rsp+\I4], r11

    mov r11,[rsp+\I5]
    mulx rcx, rax, [rip + p + 40]
    adcx r11, rcx
    adox r11, rax
    mov [rsp+\I5], r11

    mov r11,[rsp+\I6]
    mulx rbx, rax, [rip + p + 48]
    adcx r11, rbx
    adox r11, rax
    mov [rsp+\I6], r11

    mov r11,[rsp+\I7]
    mulx rcx, rax, [rip + p + 56]
    adcx r11, rcx
    adox r11, rax
    mov [rsp+\I7], r11

    mov r11,[rsp+\I8]
    mulx rbx, rax, [rip + p + 64]
    adcx r11, rbx
    adox r11, rax
    mov [rsp+\I8], r11

    mov r11,[rsp+\I9]
    mulx rcx, rax, [rip + p + 72]
    adcx r11, rcx
    adox r11, rax
    mov [rsp+\I9], r11

    mov r11,[rsp+\I10]
    mulx rbx, rax, [rip + p + 80]
    adcx r11, rbx
    adox r11, rax
    mov [rsp+\I10], r11

    mov r11,[rsp+\I11]
    mulx rcx, rax, [rip + p + 88]
    adcx r11, rcx
    adox r11, rax
    mov [rsp+\I11], r11

    mov r11,[rsp+\I12]
    mulx rbx, rax, [rip + p + 96]
    adcx r11, rbx
    adox r11, rax
    mov [rsp+\I12], r11

    mov r11,[rsp+\I13]
    mulx rcx, rax, [rip + p + 104]
    adcx r11, rcx
    adox r11, rax
    mov [rsp+\I13], r11

    mov r11,[rsp+\I14]
    mulx rbx, rax, [rip + p + 112]
    adcx r11, rbx
    adox r11, rax
    mov [rsp+\I14], r11

    mov r11,[rsp+\I15]
    mulx rcx, rax, [rip + p + 120]
    adcx r11, rcx
    adox r11, rax
    mov [rsp+\I15], r11

    mov r11,[rsp+\I16]
    mulx rbx, rax, [rip + p + 128]
    adcx r11, rbx
    adox r11, rax
    mov [rsp+\I16], r11

    mov r11,[rsp+\I17]
    mulx rcx, rax, [rip + p + 136]
    adcx r11, rcx
    adox r11, rax
    mov [rsp+\I17], r11

    mov r11,[rsp+\I18]
    mulx rbx, rax, [rip + p + 144]
    adcx r11, rbx
    adox r11, rax
    mov [rsp+\I18], r11

    mov r11,[rsp+\I19]
    mulx rcx, rax, [rip + p + 152]
    adcx r11, rcx
    adox r11, rax
    mov [rsp+\I19], r11

    mov r11,[rsp+\I20]
    mulx rbx, rax, [rip + p + 160]
    adcx r11, rbx
    adox r11, rax
    mov [rsp+\I20], r11

    mov r11,[rsp+\I21]
    mulx rcx, rax, [rip + p + 168]
    adcx r11, rcx
    adox r11, rax
    mov [rsp+\I21], r11

    mov r11,[rsp+\I22]
    mulx rbx, rax, [rip + p + 176]
    adcx r11, rbx
    adox r11, rax
    mov [rsp+\I22], r11

    mov r11,[rsp+\I23]
    mulx rcx, rax, [rip + p + 184]
    adcx r11, rcx
    adox r11, rax
    mov [rsp+\I23], r11

    mov r11,[rsp+\I24]
    mulx rbx, rax, [rip + p + 192]
    adcx r11, rbx
    adox r11, rax
    mov [rsp+\I24], r11

    mov rdx, [rdi + 8*k]

    xor rax, rax

    mov r11,[rsp+\I0]
    mulx rbx, rax, [rsi +  0]
    adox r11, rax
    mov [rsp+\I0], r11

    mov r11,[rsp+\I1]
    mulx rcx, rax, [rsi + 8]
    adcx r11, rcx
    adox r11, rax
    mov [rsp+\I1], r11

    mov r11,[rsp+\I2]
    mulx rbx, rax, [rsi + 16]
    adcx r11, rbx
    adox r11, rax
    mov [rsp+\I2], r11

    mov r11,[rsp+\I3]
    mulx rcx, rax, [rsi + 24]
    adcx r11, rcx
    adox r11, rax
    mov [rsp+\I3], r11

    mov r11,[rsp+\I4]
    mulx rbx, rax, [rsi + 32]
    adcx r11, rbx
    adox r11, rax
    mov [rsp+\I4], r11

    mov r11,[rsp+\I5]
    mulx rcx, rax, [rsi + 40]
    adcx r11, rcx
    adox r11, rax
    mov [rsp+\I5], r11

    mov r11,[rsp+\I6]
    mulx rbx, rax, [rsi + 48]
    adcx r11, rbx
    adox r11, rax
    mov [rsp+\I6], r11

    mov r11,[rsp+\I7]
    mulx rcx, rax, [rsi + 56]
    adcx r11, rcx
    adox r11, rax
    mov [rsp+\I7], r11

    mov r11,[rsp+\I8]
    mulx rbx, rax, [rsi + 64]
    adcx r11, rbx
    adox r11, rax
    mov [rsp+\I8], r11

    mov r11,[rsp+\I9]
    mulx rcx, rax, [rsi + 72]
    adcx r11, rcx
    adox r11, rax
    mov [rsp+\I9], r11

    mov r11,[rsp+\I10]
    mulx rbx, rax, [rsi + 80]
    adcx r11, rbx
    adox r11, rax
    mov [rsp+\I10], r11

    mov r11,[rsp+\I11]
    mulx rcx, rax, [rsi + 88]
    adcx r11, rcx
    adox r11, rax
    mov [rsp+\I11], r11

    mov r11,[rsp+\I12]
    mulx rbx, rax, [rsi + 96]
    adcx r11, rbx
    adox r11, rax
    mov [rsp+\I12], r11

    mov r11,[rsp+\I13]
    mulx rcx, rax, [rsi + 104]
    adcx r11, rcx
    adox r11, rax
    mov [rsp+\I13], r11

    mov r11,[rsp+\I14]
    mulx rbx, rax, [rsi + 112]
    adcx r11, rbx
    adox r11, rax
    mov [rsp+\I14], r11

    mov r11,[rsp+\I15]
    mulx rcx, rax, [rsi + 120]
    adcx r11, rcx
    adox r11, rax
    mov [rsp+\I15], r11

    mov r11,[rsp+\I16]
    mulx rbx, rax, [rsi + 128]
    adcx r11, rbx
    adox r11, rax
    mov [rsp+\I16], r11

    mov r11,[rsp+\I17]
    mulx rcx, rax, [rsi + 136]
    adcx r11, rcx
    adox r11, rax
    mov [rsp+\I17], r11

    mov r11,[rsp+\I18]
    mulx rbx, rax, [rsi + 144]
    adcx r11, rbx
    adox r11, rax
    mov [rsp+\I18], r11

    mov r11,[rsp+\I19]
    mulx rcx, rax, [rsi + 152]
    adcx r11, rcx
    adox r11, rax
    mov [rsp+\I19], r11

    mov r11,[rsp+\I20]
    mulx rbx, rax, [rsi + 160]
    adcx r11, rbx
    adox r11, rax
    mov [rsp+\I20], r11

    mov r11,[rsp+\I21]
    mulx rcx, rax, [rsi + 168]
    adcx r11, rcx
    adox r11, rax
    mov [rsp+\I21], r11

    mov r11,[rsp+\I22]
    mulx rbx, rax, [rsi + 176]
    adcx r11, rbx
    adox r11, rax
    mov [rsp+\I22], r11

    mov r11,[rsp+\I23]
    mulx rcx, rax, [rsi + 184]
    adcx r11, rcx
    adox r11, rax
    mov [rsp+\I23], r11

    mov r11,[rsp+\I24]
    mulx rbx, rax, [rsi + 192]
    adcx r11, rbx
    adox r11, rax
    mov [rsp+\I24], r11

.endm

MULSTEP 0, 0, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 136, 144, 152, 160, 168, 176, 184, 192
MULSTEP 1, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 136, 144, 152, 160, 168, 176, 184, 192, 0
MULSTEP 2, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 136, 144, 152, 160, 168, 176, 184, 192, 0, 8
MULSTEP 3, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 136, 144, 152, 160, 168, 176, 184, 192, 0, 8, 16
MULSTEP 4, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 136, 144, 152, 160, 168, 176, 184, 192, 0, 8, 16, 24
MULSTEP 5, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 136, 144, 152, 160, 168, 176, 184, 192, 0, 8, 16, 24, 32
MULSTEP 6, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 136, 144, 152, 160, 168, 176, 184, 192, 0, 8, 16, 24, 32, 40
MULSTEP 7, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 136, 144, 152, 160, 168, 176, 184, 192, 0, 8, 16, 24, 32, 40, 48
MULSTEP 8, 64, 72, 80, 88, 96, 104, 112, 120, 128, 136, 144, 152, 160, 168, 176, 184, 192, 0, 8, 16, 24, 32, 40, 48, 56
MULSTEP 9, 72, 80, 88, 96, 104, 112, 120, 128, 136, 144, 152, 160, 168, 176, 184, 192, 0, 8, 16, 24, 32, 40, 48, 56, 64
MULSTEP 10, 80, 88, 96, 104, 112, 120, 128, 136, 144, 152, 160, 168, 176, 184, 192, 0, 8, 16, 24, 32, 40, 48, 56, 64, 72
MULSTEP 11, 88, 96, 104, 112, 120, 128, 136, 144, 152, 160, 168, 176, 184, 192, 0, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80
MULSTEP 12, 96, 104, 112, 120, 128, 136, 144, 152, 160, 168, 176, 184, 192, 0, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88
MULSTEP 13, 104, 112, 120, 128, 136, 144, 152, 160, 168, 176, 184, 192, 0, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96
MULSTEP 14, 112, 120, 128, 136, 144, 152, 160, 168, 176, 184, 192, 0, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104
MULSTEP 15, 120, 128, 136, 144, 152, 160, 168, 176, 184, 192, 0, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112
MULSTEP 16, 128, 136, 144, 152, 160, 168, 176, 184, 192, 0, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120
MULSTEP 17, 136, 144, 152, 160, 168, 176, 184, 192, 0, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128
MULSTEP 18, 144, 152, 160, 168, 176, 184, 192, 0, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 136
MULSTEP 19, 152, 160, 168, 176, 184, 192, 0, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 136, 144
MULSTEP 20, 160, 168, 176, 184, 192, 0, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 136, 144, 152
MULSTEP 21, 168, 176, 184, 192, 0, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 136, 144, 152, 160
MULSTEP 22, 176, 184, 192, 0, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 136, 144, 152, 160, 168
MULSTEP 23, 184, 192, 0, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 136, 144, 152, 160, 168, 176
MULSTEP 24, 192, 0, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 136, 144, 152, 160, 168, 176, 184

    mov rdi,[rsp+208]
    mov r11,[rsp+0]
    mov [rdi+0], r11
    mov r11,[rsp+8]
    mov [rdi+8], r11
    mov r11,[rsp+16]
    mov [rdi+16], r11
    mov r11,[rsp+24]
    mov [rdi+24], r11
    mov r11,[rsp+32]
    mov [rdi+32], r11
    mov r11,[rsp+40]
    mov [rdi+40], r11
    mov r11,[rsp+48]
    mov [rdi+48], r11
    mov r11,[rsp+56]
    mov [rdi+56], r11
    mov r11,[rsp+64]
    mov [rdi+64], r11
    mov r11,[rsp+72]
    mov [rdi+72], r11
    mov r11,[rsp+80]
    mov [rdi+80], r11
    mov r11,[rsp+88]
    mov [rdi+88], r11
    mov r11,[rsp+96]
    mov [rdi+96], r11
    mov r11,[rsp+104]
    mov [rdi+104], r11
    mov r11,[rsp+112]
    mov [rdi+112], r11
    mov r11,[rsp+120]
    mov [rdi+120], r11
    mov r11,[rsp+128]
    mov [rdi+128], r11
    mov r11,[rsp+136]
    mov [rdi+136], r11
    mov r11,[rsp+144]
    mov [rdi+144], r11
    mov r11,[rsp+152]
    mov [rdi+152], r11
    mov r11,[rsp+160]
    mov [rdi+160], r11
    mov r11,[rsp+168]
    mov [rdi+168], r11
    mov r11,[rsp+176]
    mov [rdi+176], r11
    mov r11,[rsp+184]
    mov [rdi+184], r11
    mov r11,[rsp+192]
    mov [rdi+192], r11
    add rsp,216
    pop rbx
    pop rbp
    jmp .reduce_once

.global fp_mul2
fp_mul2:
    mov rdx, rsi
    mov rsi, rdi
    jmp fp_mul3

.global fp_enc
fp_enc:
    lea rdx, [rip + r_squared_mod_p]
    jmp fp_mul3

.global fp_dec
fp_dec:
    lea rdx, [rip + uint_custom_1]
    jmp fp_mul3
